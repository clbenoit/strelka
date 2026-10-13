#!/usr/bin/env python3
"""
pyflow threading smoke test (change strelka-pyflow-py3-threading-fix).

Regression test for the CPython private-name collision in the vendored
pyflow 1.1.20: StoppableThread.__init__ set `self._stop =
threading.Event()`, which shadows the private Thread._stop() method that
CPython 3.8-3.12 defines (invoked by _wait_for_tstate_lock(), reached
from is_alive() on 3.12 and from join() on 3.8-3.11) ->
`TypeError: 'Event' object is not callable` on the first harvested task
thread (observed on CPython 3.12.9, somatic/germline demos, 2026-10-13,
medfront109).

The test drives the exact CPython entry points (start -> complete ->
is_alive()/join()) on a StoppableThread, so it fails on pre-fix trees
under CPython 3.8-3.12 and passes on all supported interpreters
post-fix. It also exercises the pyflow stop API (stop/stopAll/stopped)
whose internals were renamed (self._stop -> self._stopEvent).

Usage:
    python3 pyflow_threading_smoke.py <strelka-install-tree>
The install tree must provide lib/python/pyflow/pyflow.py (any python3).
"""

import os.path
import sys
import time


def main(installDir):
    pyflowDir = os.path.join(installDir, "lib", "python", "pyflow")
    sys.path.append(pyflowDir)
    import pyflow
    if not os.path.samefile(pyflow.__file__, os.path.join(pyflowDir, "pyflow.py")):
        sys.exit("error: imported pyflow is not the installed one (%s)"
                 % pyflow.__file__)
    StoppableThread = pyflow.StoppableThread

    class Worker(StoppableThread):
        def __init__(self):
            StoppableThread.__init__(self)
            self.ran = False

        def run(self):
            time.sleep(0.1)
            self.ran = True

    # 1. lifecycle through the CPython entry points that crashed pre-fix:
    #    start -> (thread completes) -> join() -> is_alive()
    t = Worker()
    t.start()
    t.join(10)
    assert t.ran, "worker did not run"
    assert not t.is_alive(), "thread still alive after join"

    # 2. pyflow stop API (renamed internals must keep working)
    t2 = Worker()
    t2.stop()
    assert t2.stopped(), "stopped() not set after stop()"
    StoppableThread.stopAll()
    t3 = Worker()
    assert t3.stopped(), "stopped() not set after stopAll()"

    # 3. the demo pattern that crashed on 3.12.9: several short-lived
    #    threads harvested in a row (is_alive on each completed thread)
    for _ in range(10):
        w = Worker()
        w.start()
        w.join(10)
        assert w.ran
        assert not w.is_alive()

    print("pyflow threading smoke: OK (python %s)" % sys.version.split()[0])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
