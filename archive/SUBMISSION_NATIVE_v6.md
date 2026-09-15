# Measured lattice-gas vortex fluctuations and a testable fluid-model bridge

Team Merlin Digital | Airbus | Native physics revision 6 | 10 September 2026

## Delivered result

IBM Marrakesh job d9solsgpdb6s73e6a6q0 used 42 qubits, or 21 coupled pairs,
with a Taylor-Green-shaped parity preparation and an interacting lattice-gas
step. Retrieved counts support the P4 connected fluctuation observable X(q,k).
The package retains its mean-field/QTT result as a separate solver result.
This is measured quantum-medium dynamics; it is not an incompressible
Navier-Stokes solution or a demonstrated improvement over every classical solver.

## Operator and observable, kept exact

The local coupled pair uses kappa D-sum X with structural kappa. Inter-pair XX+YY
exchange supplies transport. The ordered RX/ZZ/exchange gates define a discrete
Floquet step. They are not all commuting, so exp(-iHt) is only a continuous-time
interpretation after a step-size convergence check. The revised target is the
actual emitted step and its measured covariance.

P4 is the Fourier transform in space of equal-time connected pair covariance.
Varying k gives a time-dependent structure factor. A Fourier transform over k
would give the temporal content of that statistic; it is not automatically the
two-time dynamic structure factor of a fluid. This distinction fixes the measured
quantity before an advantage is claimed.

## New audit and consequence

The retrieved P4 depth grids are nonuniform. The shared audit rejects an ordinary
FFT on those grids. Interpolation cannot manufacture new measured time points;
any nonuniform fit needs its estimator, conditioning and uncertainty reported.
The earlier retrieved P4 result remains useful as covariance evidence, but a
uniform device spectrum is a separate experiment.

A mean Taylor-Green mode has an inexpensive analytic/classical representation.
Its inability to reproduce quantum revivals compares different physical targets
and is not a same-task computational advantage. A sister experiment's operator
growth also cannot serve as a lower bound on this circuit's observable.

## Native route to a stronger application result

Use paired initial perturbations and connected covariance readout under the same
Floquet step. Test no-exchange, decoupled-pair and phase-randomized controls, then
increase lattice width at fixed local depth. Independent disconnected tiles are
parallel samples, not increased interacting system size. A two-time response may
be introduced as a new protocol with an explicit probe, rather than relabelling P4.

The supplied prepare_application_series.py emits 68 logical circuits: two prepared
modes, each with coupled and decoupled arms at every k=0..16. The campaign is
557,056 shots per run at 8,192 shots each; the largest circuit has 1,952 logical
CX before routing. The cost is explicit and the long-depth points are not yet
device-qualified. grade_application_series.py extracts the same off-diagonal
connected covariance with finite-shot bias correction; uncertainties and
independent repeats remain necessary for significance claims.

The first follow-up samples every integer k=0..16 for the same prepared modes,
with matched shots, controls and independent replication. Freeze spatial q,
temporal estimator, shot covariance and error bands before the run. Native
real-time sampling has no continuation step, but direct classical simulation of
that same discrete circuit remains an eligible comparator.

In parallel, validate the connection to the Airbus fluid objective: compare the
model's transport and damping against the stated PDE reference on held-out
parameters, with errors and any scaling map explicit. Only after that bridge
passes can a hardware response be credited as a fluid-screening benefit.

## Acceptance standard

The physics program fixes the sector, boundary treatment, initial state, operator,
probe, observable and unit map before compilation. Structural kappa is fixed by the
framework; application geometry and scientific interpretation are explicit inputs.
Fixed constants remove refitting of those constants, not the need to validate a
new application or operating condition.

The primary comparison uses the same instance and observable, the same accuracy
target, and total preparation, compilation, execution, readout and post-processing
cost. Exact small-instance calculations anchor the implementation. Direct classical
real-time methods, converged tensor networks and application-specific algorithms
remain eligible comparators. Analytic-continuation failure alone does not exclude
them. Hardware superiority is reported only after its confidence bounds separate
from those of the best matched comparator, with controls and independent replication.

## Reproducibility

This report supersedes earlier narrative claims; original reports are retained in
legacy_reports/. Existing experimental receipts keep their original identity. New
classical model checks are labelled as such and never become cloud results. The
package includes NATIVE_EXPERIMENT.json, EVIDENCE_STATUS.json and the shared
FRAMEWORK_AND_ADVANTAGE.md. Run the shared validate_upgrades.py from the parent
directory to verify integrity and the new mathematical checks. No new hardware
measurement is claimed by this revision.
