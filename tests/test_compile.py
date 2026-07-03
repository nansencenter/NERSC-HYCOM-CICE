"""
Phase 2 tests — compile hycom_cice for TP0 expt_01.0.

These tests require either:
  • Native Linux with conda-forge ESMF in the environment, or
  • macOS with Podman and the ``hycom-esmf-spike`` image built.

They are skipped automatically when neither is available (e.g. a developer
machine with none of the build toolchain set up).
"""


def test_compile_succeeds(compiled_hycom):
    """compile_model.sh must exit 0 for TP0 expt_01.0 with gfortran+OpenMPI+ESMF."""
    assert compiled_hycom["success"], (
        "compile_model.sh failed.\n"
        "--- stdout ---\n" + compiled_hycom["stdout"][-3000:] + "\n"
        "--- stderr ---\n" + compiled_hycom["stderr"][-3000:]
    )


def test_compile_produces_executable(compiled_hycom):
    """The build log must mention creation of the hycom_cice executable."""
    assert compiled_hycom["success"], "Build failed — see test_compile_succeeds for details"
    combined = compiled_hycom["stdout"] + compiled_hycom["stderr"]
    assert "hycom_cice" in combined, \
        "Expected 'hycom_cice' to appear in build output"
