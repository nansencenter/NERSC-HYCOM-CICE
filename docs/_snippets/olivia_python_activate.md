Prepend the `bin` directory of your containerised environment to `PATH`. Replace
`<install_dir>` with the path you used during installation:

```bash
export PATH="<install_dir>/bin:${PATH}"
```

This line should be present in `~/.bashrc` and in any job submission script that uses
the model's Python tools.
