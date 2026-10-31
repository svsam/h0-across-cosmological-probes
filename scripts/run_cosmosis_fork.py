#!/usr/bin/env python3
"""Run the existing CosmoSIS installation with POSIX fork for local SMP."""

import multiprocessing as multiprocessing

multiprocessing.set_start_method("fork", force=True)

import cosmosis.main
from cosmosis.runtime import process_pool

# This local CosmoSIS version's multiprocessing pool has placeholder
# collectives that return an unset attribute after sampling. The sampler's
# shared-memory pool is used only on the master process for these calls.
process_pool.Pool.bcast = lambda self, data: data
process_pool.Pool.gather = lambda self, data: data
process_pool.Pool.allreduce = lambda self, data: data


def persistent_map(self, function, args):
    # Create workers after the sampler has initialized its global pipeline.
    if not hasattr(self, "_workers"):
        self._workers = multiprocessing.Pool(self.size)
    return self._workers.map(function, args)


def close_workers(self, *args):
    if hasattr(self, "_workers"):
        self._workers.close()
        self._workers.join()


process_pool.Pool.map = persistent_map
process_pool.Pool.close = close_workers
process_pool.Pool.__exit__ = close_workers


if __name__ == "__main__":
    raise SystemExit(cosmosis.main.main())
