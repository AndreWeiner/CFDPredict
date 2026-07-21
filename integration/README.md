# integration/ — HPC / Slurm / online integration (TU Dresden)

Maintained by **TU Dresden**. Slurm wrappers, HPC submission, and online ROM
integration around the DHCAE workflows — which stay **unmodified**.

`slurm/run_a2.sbatch.example` is a reference batch script contributed by DHCAE from
a 2-node Slurm sandbox test, where the A2 workflow ran end-to-end unchanged
(cross-node, identical objectives to the 128-core reference run).

> **Ownership boundary:** this directory belongs to TU Dresden. The DHCAE workflow
> code under `workflows/`, `tools/`, `examples/`, `docs/` is wrapped, not edited.
