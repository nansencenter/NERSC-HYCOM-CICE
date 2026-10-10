"""
Compilation tests — build hycom_cice for TP0 expt_01.0.

Requires either native Linux with conda-forge ESMF, or macOS with Podman
and the ``hycom-esmf-spike`` image built. Skipped automatically when neither
is available.
"""

from typing import Any


def test_compile_succeeds(compiled_hycom: dict[str, Any]) -> None:
    """compile_model.sh must exit 0 for TP0 expt_01.0 with gfortran+OpenMPI+ESMF."""
    assert compiled_hycom["success"], (
        "compile_model.sh failed.\n"
        "--- stdout ---\n" + compiled_hycom["stdout"][-3000:] + "\n"
        "--- stderr ---\n" + compiled_hycom["stderr"][-3000:]
    )


def test_compile_produces_executable(compiled_hycom: dict[str, Any]) -> None:
    """The build log must mention creation of the hycom_cice executable."""
    assert compiled_hycom["success"], "Build failed — see test_compile_succeeds for details"
    combined = compiled_hycom["stdout"] + compiled_hycom["stderr"]
    assert "hycom_cice" in combined, \
        "Expected 'hycom_cice' to appear in build output"
