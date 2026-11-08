"""CPU-only AlphaFold 3 container wrapper."""

import sys
import jax
import jax.nn
import jax.numpy as jnp
import tokamax


def _cpu_dot_product_attention(
    query, key, value, bias=None, mask=None, **kwargs
):
  del kwargs
  scale = 1.0 / jnp.sqrt(query.shape[-1])
  logits = jnp.einsum("...qhd,...khd->...hqk", query, key) * scale
  if bias is not None:
    logits = logits + bias
  if mask is not None:
    logits = jnp.where(mask, logits, -1e9)
  weights = jax.nn.softmax(logits, axis=-1)
  return jnp.einsum("...hqk,...khd->...qhd", weights, value)


tokamax.dot_product_attention = _cpu_dot_product_attention
_orig = jax.local_devices
jax.local_devices = (
    lambda process_index=None, backend=None, host_id=None: _orig(
        process_index=process_index, backend="cpu", host_id=host_id
    )
)
with open("/app/alphafold/run_alphafold.py", "r", encoding="utf-8") as f:
  code = f.read()
code = code.replace(
    "gpu_devices[_GPU_DEVICE.value].compute_capability", "'8.0'"
)
sys.argv = ["/app/alphafold/run_alphafold.py"] + sys.argv[1:]
sys.path.insert(0, "/app/alphafold")
exec(
    compile(code, "/app/alphafold/run_alphafold.py", "exec"),
    {"__name__": "__main__", "__file__": "/app/alphafold/run_alphafold.py"},
)
