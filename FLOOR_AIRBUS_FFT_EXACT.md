# Floor: the v8 time-to-solution crossover against the exact FFT solve (laptop, 2026-09-14 15:50)

**What v8 claims (§1, §4, App. C).** Fault-tolerant Schrödingerisation of A = −U_c D_x + νL overtakes
classical explicit finite volume at N ≈ 65,536 per axis at every Re, 61–71 logical qubits; business
scenario "32 days against 51 days".

**The cheapest adversary.** A is circulant on the periodic cell grid, so exp(AT)u₀ is one 2D FFT, one
diagonal multiply, one inverse FFT. No time stepping, O(N² log N), and it returns the very vector the
Schrödingerised readout samples. The package already knows this: `schrod_tgv.py` computes its
"classical FV" validation column exactly this way (line ≈ 78, comment "exact classical solve of the same
operator") and simulates H_S by diagonalising it in the DFT basis. The crossover in §4 is measured against
a different, slower classical algorithm (explicit stepping, N²·steps flops) that nobody would use here.

**Measured (this laptop, `floor_airbus_fft_exact.py` → `floor_airbus_fft_exact.json`, T = 0.5).**

| Re | N | FFT exact solve | L2 vs analytic (package convention) | v8 quantum time at 1 MHz | v8 explicit-FV model at 10 GFLOP/s |
|---:|---:|---:|---:|---:|---:|
| 100 | 128 | 3 ms | 1.34×10⁻⁴ (= v8 table) | 214 s | 6 ms |
| 100 | 4,096 | 7.1 s | 1.3×10⁻⁷ | 6.2×10⁵ s | 6.7×10³ s |
| 10⁴ | 128 | 2 ms | 1.42×10⁻⁴ (= v8 table) | 2.7 s | 4 ms |
| 10⁴ | 4,096 | 6.2 s | 1.4×10⁻⁷ | 6.1×10³ s | 131 s |

Extrapolated with the measured N² log N fit to the claimed crossover N = 65,536: FFT ≈ 2,250 s (37 min,
69 GB of complex field, or seconds and 128 bytes in the spectral representation since the flow is one
mode); v8's quantum time 2.75×10⁶ s (Re 10⁴) to 2.77×10⁹ s (Re 10): **1.2×10³ to 1.2×10⁶ times slower
than the classical exact solve** at the rung v8 calls the crossing, and slower at every smaller rung too.

**Verdict.** The v8 time-to-solution advantage does not exist against the correct classical comparator,
and the correct comparator is inside the package. The memory claim was already conceded against the
128-byte spectral floor in §8(ii). What v8 leaves standing on advantage: nothing measured or modelled.
The v7 tensor-train spine had at least a memory ratio against dense storage (also floored, 18–105× behind
spectral). Grade on the advantage criterion: v8 is below v7.

**What a crossing on this track would need.** A linear operator that is *not* diagonalisable by FFT
(non-constant coefficients: linearised stability about a non-uniform base flow, variable-viscosity or
mapped-grid advection–diffusion), where the classical exact route is Krylov at O(N²·steps) and the
Schrödingerised ‖H‖T·polylog N scaling is a genuine polynomial gap; or a parametric ensemble (Re, initial
condition, or geometry registers in superposition) where classical cost multiplies by the ensemble size and
the quantum register grows by its log; or the nonlinear multi-mode regime through Carleman, inside the
O(Re^{3D/8}) bound v8 already cites. The canonical vortex itself gives no wall, and a submission that says
so in one line and puts the crossing on a named non-circulant Airbus operator is the honest A-path.
Nothing here needs the workstation; the decisive computation took seven seconds.
