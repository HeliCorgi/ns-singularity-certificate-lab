# Hou early-time critical-`L^3` analysis

## Scope

This report was computed from the checkpoint files in the last repository ZIP available in this conversation (`65x128`, `129x256`, `193x384`, five times through `T1 = 0.002191729`).  It is not a run of the newer 711-test GitHub branch.  The same post-processor must be rerun on the newer `257x512` checkpoint and on every future whole-space candidate.

All quantities below are floating-point observations on the represented finite periodic cylinder.  They are not `L^3(R^3)` bounds and are not a proof of regularity or singularity.

## Result

| grid | `L3(T1)/L3(0)` | `max|u|(T1)/max|u|(0)` | `Lr: 0 -> T1` | `Lz: 0 -> T1` | `Q(T1)/Q(0)` | effective shells `0 -> T1` | outer radial fraction at `T1` |
|---|---:|---:|---:|---:|---:|---:|---:|
| 65x128 | 0.901469 | 1.308380 | 0.208514 -> 0.260737 | 0.139905 -> 0.169170 | 4.234727 | 1.914142 -> 1.771682 | 1.170e-07 |
| 129x256 | 0.908517 | 1.876747 | 0.208514 -> 0.258955 | 0.139905 -> 0.166939 | 12.165162 | 1.916061 -> 1.803189 | 5.328e-08 |
| 193x384 | 0.910394 | 2.071700 | 0.208514 -> 0.258422 | 0.139905 -> 0.166371 | 16.240893 | 1.915033 -> 1.813706 | 3.401e-08 |

At the finest available grid, the represented-domain critical norm changed as

```text
L3: 118.565484636 -> 107.941332825
max |u|: 327.785671 -> 679.073580
Lr: 0.20851441 -> 0.25842155
Lz: 0.13990545 -> 0.16637138
```

The physical velocity maximum increased, but both critical-density RMS widths **expanded**, rather than contracted.  The global represented-domain `L3` norm decreased by about `8.96%`.  Therefore the large increase in the RMS diagnostic `Q = A^3 Lr^2 Lz` is not evidence of an anisotropic Type-II concentrating profile; it is dominated by the growing pointwise amplitude while the global critical density remains broad.

The effective dyadic shell count decreased slightly, from about `1.915` to `1.814`, and the outer radial critical-mass fraction remained about `3.4e-8` on the finest grid.  In this early periodic-cylinder window there is no observed growth of global critical mass, no observed proliferation of occupied scales, and no observed radial outer-tail accumulation.

The first moving shell also lost absolute critical mass on all three grids:

| grid | first-shell mass at `t=0` | first-shell mass at `T1` | ratio |
|---|---:|---:|---:|
| 65x128 | 315817 | 232438 | 0.735989 |
| 129x256 | 315450 | 248107 | 0.786516 |
| 193x384 | 315083 | 253534 | 0.804659 |

Because the shell radius is defined from the RMS widths and those widths expand, this last observation is only a diagnostic.  It nevertheless gives no early support for critical-mass concentration.

## Scientific conclusion

The early Hou run remains a resolved-enough observation of strong vorticity amplification and changing local geometry.  In the critical topology required to evade the endpoint regularity obstruction, however, this time window is negative:

> The available early snapshots show decreasing finite-cylinder `L3`, expanding critical-density widths, nearly constant or decreasing shell complexity, and negligible radial outer mass.  They do not yet exhibit Type-II critical growth, anisotropic critical concentration, a multiscale cascade, or an outer-tail mechanism.

This does not prove that the later Hou evolution is regular.  It shows that the current early-time data have not crossed the necessary critical-norm gate.

## Required next computation

1. Run the same diagnostics on the latest `257x512` snapshots.
2. Integrate the diagnostics into full-step streaming output rather than checkpoint-only post-processing.
3. Add nonperiodic-`z` compactly supported data and a free-space elliptic solver.
4. Repeat under independent `Rmax` and `Zmax` enlargement.
5. Track fixed physical shells and moving rescaled shells simultaneously.
6. Reject a whole-space candidate if a validated global `L3` bound remains uniform.
7. Promote a candidate only after critical growth survives space, time, domain and integrator refinement.
