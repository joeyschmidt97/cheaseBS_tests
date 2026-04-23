# Iterative CHEASE Workflow

- Boundary is fixed to the original EQDSK boundary.
- Target Ip is fixed to the original EQDSK value.
- Pressure uses `p_total = p_thermal(new) + p_fast(original floored)`.
- Bootstrap and driven current are combined on the `j_parallel` footing, then converted to total `I*` before each replay solve.
- Primary replay coordinate is `rhop`.
- Bootstrap is extracted diagnostically from each solved equilibrium via `ogyropsi.dat` and `JBSBAV/B0`.
- Internal CHEASE bootstrap feedback is disabled during the equilibrium solve.
- The next iteration uses `j_total_input = A * j_driven_template + j_bootstrap(previous)`.
- Final iteration: `0`.
- Final Ip: `743173.038200 A` vs target `753908.302000 A`.
- q comparison is reported at `rhop=0.985` rather than the singular last grid point.
- q(`0.985`): target `10.080731196`, final `10.133286233`.
- Final amplitude: `1.00000000`.
- Converged: `False`.
