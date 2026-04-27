use std::ffi::c_void;
use std::mem::{size_of, size_of_val};
use std::ptr::NonNull;
use std::sync::Arc;

use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::{SendableRecordBatchStream, TaskContext};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use futures::StreamExt;
use futures::stream::BoxStream;
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSString;
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandEncoder, MTLCommandQueue, MTLComputeCommandEncoder,
    MTLComputePipelineState, MTLCreateSystemDefaultDevice, MTLDevice, MTLLibrary,
    MTLResourceOptions, MTLSize,
};

use crate::count_overlaps_rank::{
    CountOverlapsRankIndex, MISSING_CONTIG_ID, append_count_column, project_batch,
};
use crate::filter_op::FilterOp;

#[link(name = "CoreGraphics", kind = "framework")]
unsafe extern "C" {}

const KERNEL_SOURCE: &str = r#"
#include <metal_stdlib>
using namespace metal;

struct ContigRange {
    uint start_offset;
    uint end_offset;
    uint len;
    uint pad;
};

static inline uint upper_bound_i64(device const long* values, uint offset, uint len, long needle) {
    uint lo = 0;
    uint hi = len;
    while (lo < hi) {
        uint mid = lo + ((hi - lo) >> 1);
        if (values[offset + mid] <= needle) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

static inline uint lower_bound_i64(device const long* values, uint offset, uint len, long needle) {
    uint lo = 0;
    uint hi = len;
    while (lo < hi) {
        uint mid = lo + ((hi - lo) >> 1);
        if (values[offset + mid] < needle) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

kernel void count_overlaps_rank(
    device const ContigRange* contigs [[buffer(0)]],
    device const long* starts [[buffer(1)]],
    device const long* ends [[buffer(2)]],
    device const uint* query_contigs [[buffer(3)]],
    device const long* query_starts [[buffer(4)]],
    device const long* query_ends [[buffer(5)]],
    device long* output [[buffer(6)]],
    constant uint& n_queries [[buffer(7)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= n_queries) {
        return;
    }

    uint contig_id = query_contigs[gid];
    long query_start = query_starts[gid];
    long query_end = query_ends[gid];
    if (contig_id == 0xffffffffu || query_start > query_end) {
        output[gid] = 0;
        return;
    }

    ContigRange range = contigs[contig_id];
    uint started = upper_bound_i64(starts, range.start_offset, range.len, query_end);
    uint ended_prior = lower_bound_i64(ends, range.end_offset, range.len, query_start);
    output[gid] = long(started - ended_prior);
}
"#;

type MetalDevice = ProtocolObject<dyn MTLDevice>;
type MetalCommandQueue = ProtocolObject<dyn MTLCommandQueue>;
type MetalComputePipelineState = ProtocolObject<dyn MTLComputePipelineState>;
type MetalBuffer = ProtocolObject<dyn MTLBuffer>;

pub struct AppleGpuCountOverlapsBackend {
    index: CountOverlapsRankIndex,
    runtime: MetalRuntime,
    contig_buffer: SharedMetalBuffer,
    starts_buffer: SharedMetalBuffer,
    ends_buffer: SharedMetalBuffer,
}

struct MetalRuntime {
    device: Retained<MetalDevice>,
    queue: Retained<MetalCommandQueue>,
    pipeline: Retained<MetalComputePipelineState>,
}

struct SharedMetalBuffer(Retained<MetalBuffer>);

// MTLBuffer is safe to share here because these buffers are immutable after
// initialization, command buffers retain referenced resources, and batch-local
// writable buffers are never stored in this shared wrapper.
unsafe impl Send for SharedMetalBuffer {}
unsafe impl Sync for SharedMetalBuffer {}

impl SharedMetalBuffer {
    fn as_buffer(&self) -> &MetalBuffer {
        &self.0
    }
}

impl AppleGpuCountOverlapsBackend {
    pub fn try_new(index: CountOverlapsRankIndex) -> Result<Self> {
        if index.interval_count() == 0 {
            return Err(DataFusionError::Execution(
                "Apple GPU count_overlaps backend requires a non-empty left index".to_string(),
            ));
        }

        let runtime = MetalRuntime::new()?;
        let contig_buffer = runtime.new_buffer_from_slice(index.contigs())?;
        let starts_buffer = runtime.new_buffer_from_slice(index.starts_sorted())?;
        let ends_buffer = runtime.new_buffer_from_slice(index.ends_sorted())?;

        Ok(Self {
            index,
            runtime,
            contig_buffer: SharedMetalBuffer(contig_buffer),
            starts_buffer: SharedMetalBuffer(starts_buffer),
            ends_buffer: SharedMetalBuffer(ends_buffer),
        })
    }

    fn execute_batch(
        &self,
        batch: &datafusion::arrow::array::RecordBatch,
        full_schema: SchemaRef,
        output_schema: SchemaRef,
        projection: Option<&[usize]>,
        columns: (&str, &str, &str),
        filter_op: FilterOp,
    ) -> Result<datafusion::arrow::array::RecordBatch> {
        let queries = self.index.encode_query_batch(batch, columns, filter_op)?;
        if queries.query_contigs.is_empty() {
            let batch = append_count_column(batch, full_schema, Vec::new())?;
            return project_batch(batch, output_schema, projection);
        }

        let query_contigs = self.runtime.new_buffer_from_slice(&queries.query_contigs)?;
        let query_starts = self.runtime.new_buffer_from_slice(&queries.query_starts)?;
        let query_ends = self.runtime.new_buffer_from_slice(&queries.query_ends)?;
        let output = self
            .runtime
            .new_buffer_with_len(queries.query_contigs.len() * size_of::<i64>())?;

        let command_buffer = self.runtime.queue.commandBuffer().ok_or_else(|| {
            DataFusionError::Execution("failed to create Metal command buffer".to_string())
        })?;
        let encoder = command_buffer.computeCommandEncoder().ok_or_else(|| {
            DataFusionError::Execution("failed to create Metal compute encoder".to_string())
        })?;

        encoder.setComputePipelineState(&self.runtime.pipeline);
        unsafe {
            encoder.setBuffer_offset_atIndex(Some(self.contig_buffer.as_buffer()), 0, 0);
            encoder.setBuffer_offset_atIndex(Some(self.starts_buffer.as_buffer()), 0, 1);
            encoder.setBuffer_offset_atIndex(Some(self.ends_buffer.as_buffer()), 0, 2);
            encoder.setBuffer_offset_atIndex(Some(&query_contigs), 0, 3);
            encoder.setBuffer_offset_atIndex(Some(&query_starts), 0, 4);
            encoder.setBuffer_offset_atIndex(Some(&query_ends), 0, 5);
            encoder.setBuffer_offset_atIndex(Some(&output), 0, 6);

            let n_queries = u32::try_from(queries.query_contigs.len()).map_err(|_| {
                DataFusionError::Execution(
                    "count_overlaps GPU query batch exceeds u32::MAX rows".to_string(),
                )
            })?;
            encoder.setBytes_length_atIndex(
                NonNull::from(&n_queries).cast::<c_void>(),
                size_of::<u32>(),
                7,
            );
        }

        let threads_per_group = self.runtime.threads_per_threadgroup();
        let threads_per_grid = MTLSize {
            width: queries.query_contigs.len(),
            height: 1,
            depth: 1,
        };
        encoder.dispatchThreads_threadsPerThreadgroup(threads_per_grid, threads_per_group);
        encoder.endEncoding();
        command_buffer.commit();
        command_buffer.waitUntilCompleted();

        if let Some(error) = command_buffer.error() {
            return Err(DataFusionError::Execution(format!(
                "Metal count_overlaps command failed: {error}"
            )));
        }

        let counts = unsafe {
            let ptr = output.contents().as_ptr().cast::<i64>();
            std::slice::from_raw_parts(ptr, queries.query_contigs.len()).to_vec()
        };
        let batch = append_count_column(batch, full_schema, counts)?;
        project_batch(batch, output_schema, projection)
    }
}

impl MetalRuntime {
    fn new() -> Result<Self> {
        let device = MTLCreateSystemDefaultDevice().ok_or_else(|| {
            DataFusionError::Execution("no default Metal device is available".to_string())
        })?;
        let queue = device.newCommandQueue().ok_or_else(|| {
            DataFusionError::Execution("failed to create Metal command queue".to_string())
        })?;
        let source = NSString::from_str(KERNEL_SOURCE);
        let library = device
            .newLibraryWithSource_options_error(&source, None)
            .map_err(|error| {
                DataFusionError::Execution(format!(
                    "failed to compile count_overlaps Metal library: {error}"
                ))
            })?;
        let function_name = NSString::from_str("count_overlaps_rank");
        let function = library.newFunctionWithName(&function_name).ok_or_else(|| {
            DataFusionError::Execution(
                "count_overlaps_rank Metal function was not found".to_string(),
            )
        })?;
        let pipeline = device
            .newComputePipelineStateWithFunction_error(&function)
            .map_err(|error| {
                DataFusionError::Execution(format!(
                    "failed to create count_overlaps Metal pipeline: {error}"
                ))
            })?;

        Ok(Self {
            device,
            queue,
            pipeline,
        })
    }

    fn new_buffer_from_slice<T>(&self, values: &[T]) -> Result<Retained<MetalBuffer>> {
        if values.is_empty() {
            return self.new_buffer_with_len(1);
        }
        let len = size_of_val(values);
        let ptr = NonNull::new(values.as_ptr() as *mut c_void).ok_or_else(|| {
            DataFusionError::Execution("cannot create Metal buffer from null slice".to_string())
        })?;
        unsafe {
            self.device
                .newBufferWithBytes_length_options(ptr, len, MTLResourceOptions::StorageModeShared)
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "failed to create Metal buffer with {len} bytes"
                    ))
                })
        }
    }

    fn new_buffer_with_len(&self, len: usize) -> Result<Retained<MetalBuffer>> {
        self.device
            .newBufferWithLength_options(len.max(1), MTLResourceOptions::StorageModeShared)
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "failed to allocate Metal buffer with {} bytes",
                    len.max(1)
                ))
            })
    }

    fn threads_per_threadgroup(&self) -> MTLSize {
        let width = self.runtime_thread_width().clamp(1, 256);
        MTLSize {
            width,
            height: 1,
            depth: 1,
        }
    }

    fn runtime_thread_width(&self) -> usize {
        let execution_width = self.pipeline.threadExecutionWidth();
        let max_threads = self.pipeline.maxTotalThreadsPerThreadgroup();
        let preferred = if execution_width == 0 {
            128
        } else {
            execution_width * 4
        };
        preferred.min(max_threads.max(1))
    }
}

#[allow(clippy::too_many_arguments)]
pub fn get_apple_gpu_stream(
    right_plan: Arc<dyn ExecutionPlan>,
    backend: Arc<AppleGpuCountOverlapsBackend>,
    full_schema: SchemaRef,
    output_schema: SchemaRef,
    projection: Option<Arc<Vec<usize>>>,
    columns_2: Arc<(String, String, String)>,
    filter_op: FilterOp,
    partition: usize,
    context: Arc<TaskContext>,
) -> Result<SendableRecordBatchStream> {
    let partition_stream = right_plan.execute(partition, context)?;
    let full_schema_for_closure = full_schema.clone();
    let output_schema_for_closure = output_schema.clone();
    let projection_for_closure = projection.clone();
    let iter = partition_stream.map(move |rb| {
        rb.and_then(|rb| {
            backend.execute_batch(
                &rb,
                full_schema_for_closure.clone(),
                output_schema_for_closure.clone(),
                projection_for_closure.as_deref().map(Vec::as_slice),
                (&columns_2.0, &columns_2.1, &columns_2.2),
                filter_op.clone(),
            )
        })
    });
    Ok(Box::pin(RecordBatchStreamAdapter::new(
        output_schema,
        Box::pin(iter) as BoxStream<'static, Result<datafusion::arrow::array::RecordBatch>>,
    )))
}

#[allow(dead_code)]
fn _assert_missing_contig_id_matches_kernel() {
    const _: u32 = MISSING_CONTIG_ID;
}
