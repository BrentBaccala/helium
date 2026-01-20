"""Dynamic Process Manager with Signal-Based Worker Management

This module provides a simple process manager that allows dynamic addition
and removal of worker processes via Unix signals, designed for long-running
worker tasks. Uses fork-based multiprocessing for compatibility with SageMath
and other non-standard Python environments.

Features:
    - Add workers dynamically by sending SIGUSR1
    - Remove workers dynamically by sending SIGUSR2
    - Signal-aware waiting that returns immediately on signal receipt
    - Graceful shutdown on SIGTERM/SIGINT
    - Automatic collection of terminated workers
    - Compatible with SageMath and other Python variants

Usage:
    from dynamic_pool import DynamicProcessManager
    
    def my_worker():
        while True:
            do_work()
    
    if __name__ == "__main__":
        manager = DynamicProcessManager(
            worker_func=my_worker,
            num_workers=4,
            initializer=init_func,  # optional
            initargs=()             # optional
        )
        manager.start()
        
        while manager.is_running():
            exited = manager.wait_for_events()
            # Handle worker exits

Author: Claude (Anthropic)
Date: January 2026
"""

import signal
import os
import multiprocessing as mp
import time
from threading import Event, Lock


class DynamicProcessManager:
    """Manages a pool of worker processes with dynamic resize via signals."""
    
    def __init__(self, worker_func, num_workers=None, initializer=None, initargs=()):
        """Initialize the process manager.
        
        Args:
            worker_func: The function each worker process will execute
            num_workers: Initial number of workers (default: CPU count)
            initializer: Optional function to call in each worker before worker_func
            initargs: Arguments to pass to initializer
        """
        self.worker_func = worker_func
        self.initializer = initializer
        self.initargs = initargs
        self.num_workers = num_workers or os.cpu_count() or 1
        
        self.processes = []
        self.lock = Lock()
        self.shutdown_event = Event()
        self.signal_event = Event()
        
        # Register signal handlers
        signal.signal(signal.SIGUSR1, self._handle_sigusr1)
        signal.signal(signal.SIGUSR2, self._handle_sigusr2)
        signal.signal(signal.SIGTERM, self._handle_sigterm)
        signal.signal(signal.SIGINT, self._handle_sigterm)
    
    def _worker_wrapper(self):
        """Wrapper that calls initializer then worker_func."""
        if self.initializer is not None:
            try:
                self.initializer(*self.initargs)
            except Exception as e:
                print(f"Worker {os.getpid()} initializer failed: {e}")
                return
        
        try:
            self.worker_func()
        except KeyboardInterrupt:
            pass
        except Exception as e:
            print(f"Worker {os.getpid()} failed: {e}")
    
    def _spawn_worker(self):
        """Spawn a new worker process."""
        p = mp.Process(target=self._worker_wrapper)
        p.start()
        with self.lock:
            self.processes.append(p)
        print(f"Started worker {p.pid}. Total workers: {len(self.processes)}")
        return p
    
    def _handle_sigusr1(self, signum, frame):
        """Add a new worker process when SIGUSR1 is received."""
        if self.shutdown_event.is_set():
            return
        
        self._spawn_worker()
        self.signal_event.set()
    
    def _handle_sigusr2(self, signum, frame):
        """Remove a worker process when SIGUSR2 is received."""
        if self.shutdown_event.is_set():
            return
        
        with self.lock:
            if len(self.processes) > 0:
                # Terminate the last worker
                p = self.processes[-1]
                p.terminate()
                print(f"Terminating worker {p.pid}. Remaining: {len(self.processes)-1}")
            else:
                print("No workers to remove")
        
        self.signal_event.set()
    
    def _handle_sigterm(self, signum, frame):
        """Handle termination signals."""
        print("\nReceived termination signal, shutting down...")
        self.shutdown_event.set()
        self.signal_event.set()
    
    def start(self):
        """Start the initial pool of workers."""
        print(f"Starting {self.num_workers} workers (PID: {os.getpid()})")
        print(f"Send SIGUSR1 to add workers: kill -USR1 {os.getpid()}")
        print(f"Send SIGUSR2 to remove workers: kill -USR2 {os.getpid()}")
        print(f"Send SIGTERM or Ctrl+C to shutdown")
        print()
        
        for _ in range(self.num_workers):
            self._spawn_worker()
    
    def wait_for_events(self, timeout=1.0):
        """Wait for worker exits or signals.
        
        Args:
            timeout: How long to wait before checking (seconds)
        
        Returns:
            List of PIDs that exited
        """
        # Wait for signal or timeout
        self.signal_event.wait(timeout)
        self.signal_event.clear()
        
        # Check for terminated processes
        exited = []
        with self.lock:
            alive = []
            for p in self.processes:
                if not p.is_alive():
                    p.join()  # Clean up zombie
                    exited.append(p.pid)
                    print(f"Worker {p.pid} exited. Remaining: {len(self.processes)-1}")
                else:
                    alive.append(p)
            self.processes = alive
        
        return exited
    
    def is_running(self):
        """Check if the manager should continue running."""
        if self.shutdown_event.is_set():
            return False
        
        with self.lock:
            # Continue if we have workers or haven't been told to shutdown
            return len(self.processes) > 0
    
    def shutdown(self, wait=True):
        """Shutdown all workers.
        
        Args:
            wait: If True, wait for workers to terminate gracefully
        """
        print("Shutting down all workers...")
        self.shutdown_event.set()
        
        with self.lock:
            for p in self.processes:
                if wait:
                    p.terminate()
                else:
                    p.kill()
            
            if wait:
                for p in self.processes:
                    p.join(timeout=5)
                    if p.is_alive():
                        print(f"Worker {p.pid} did not terminate, killing...")
                        p.kill()
                        p.join()
            
            self.processes.clear()
        
        print("All workers shut down")
    
    def get_worker_count(self):
        """Return the current number of workers."""
        with self.lock:
            return len(self.processes)


# Example usage
import time
import random

def long_running_worker():
    """A long-running worker function that processes items indefinitely."""
    pid = os.getpid()
    print(f"Worker {pid} started")
    
    try:
        counter = 0
        while True:
            # Simulate doing work
            time.sleep(random.uniform(0.5, 1.5))
            counter += 1
            if counter % 5 == 0:
                print(f"Worker {pid} processed {counter} items")
    except KeyboardInterrupt:
        print(f"Worker {pid} interrupted")
    except Exception as e:
        print(f"Worker {pid} error: {e}")
    finally:
        print(f"Worker {pid} shutting down")


if __name__ == "__main__":
    # Create manager starting with 2 workers
    manager = DynamicProcessManager(
        worker_func=long_running_worker,
        num_workers=2
    )
    
    manager.start()
    
    try:
        # Main loop - wait for events
        while manager.is_running():
            exited = manager.wait_for_events()
            
            # Optionally handle exited workers
            if exited:
                print(f"Workers exited: {exited}")
    
    except KeyboardInterrupt:
        print("\nCtrl+C received")
    finally:
        manager.shutdown(wait=True)