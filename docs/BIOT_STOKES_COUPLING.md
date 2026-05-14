# Biot poroelasticity and Stokes–Biot coupling

This document describes the **Biot** equation set and the closure evaluators used to couple a poroelastic porous medium to an incompressible fluid region modeled with **Navier–Stokes** (or Stokes) fields in Flujo.

The formulation follows the thick Biot layer and Stokes–Biot interface conditions in Scharf, Bukač, and Čanić, [*Splitting method for a multilayered poroelastic solid interacting with Stokes flow*](https://arxiv.org/abs/2507.10538) (arXiv:2507.10538). Flujo implements the **3D thick Biot subdomain** and **interface closures**; it does not yet include the thin poroelastic plate model from that paper.

---

## Source layout

| Path | Role |
|------|------|
| `src/eqn_sets/Flujo_EquationSet_Biot.hpp` | `EquationSet_Biot` declaration |
| `src/eqn_sets/Flujo_EquationSet_Biot_impl.hpp` | DOF registration and residual evaluator assembly |
| `src/closures/Flujo_BiotPoroelasticStress.{hpp,_impl.hpp}` | Poroelastic stress flux closure |
| `src/closures/Flujo_StokesBiotFluidTraction.{hpp,_impl.hpp}` | Fluid Cauchy traction on the interface |
| `src/closures/Flujo_StokesBiotNormalFiltration.{hpp,_impl.hpp}` | Normal velocity jump for kinematic coupling |
| `src/closures/Flujo_StokesBiotBeaversJosephSlip.{hpp,_impl.hpp}` | Beavers–Joseph–Saffman (BJS) slip residual |
| `src/closures/Flujo_ClosureModel_Factory.{hpp,_impl.hpp}` | Registers Biot and Stokes–Biot closure types |
| `src/closures/Flujo_ClosureModel_Factory_TemplateBuilder.hpp` | Panzer template builder for the closure factory |
| `src/eqn_sets/Flujo_EquationSetFactory.hpp` | Registers physics type `Biot` |

---

## Continuous model (thick Biot domain)

In the poroelastic region \(\Omega_b\), the first-order Biot system uses displacement \(\eta\), structure velocity \(\xi = \partial_t \eta\), pore pressure \(p\), and Darcy filtration velocity \(\mathbf{u}_b\):

\[
\rho_b \partial_t \xi - \nabla \cdot \sigma_b(\eta, p) + \gamma \eta = \mathbf{F}_b,
\]
\[
c_0 \partial_t p + \alpha \nabla \cdot \xi + \nabla \cdot \mathbf{u}_b = G_b,
\]
\[
\kappa^{-1} \mathbf{u}_b = -\nabla p,
\]

with poroelastic stress

\[
\sigma_b(\eta, p) = 2\mu_b D(\eta) + \lambda_b (\nabla \cdot \eta)\, I - \alpha p\, I,
\]

and \(D(\eta) = \tfrac{1}{2}(\nabla \eta + (\nabla \eta)^T)\).

On the fluid–structure interface \(\Gamma\), the reference paper couples Stokes flow to Biot through normal velocity continuity, displacement kinematics, stress balance, pressure transfer, and optional BJS slip. The closure evaluators below expose those quantities as fields for residuals or boundary strategies.

---

## Panzer registration

### Physics type

In a physics block equation set, set:

```text
Type: Biot
```

### Degrees of freedom

| Field | Default name | Default basis | Registered auxiliaries |
|-------|----------------|---------------|-------------------------|
| Solid displacement | `displacement` | `HGrad` (vector) | `GRAD_displacement`, `DXDT_displacement` (transient) |
| Pore pressure | `pore_pressure` | `HGrad` (scalar) | `GRAD_pore_pressure`, `DXDT_pore_pressure` (transient) |
| Filtration velocity | `filtration_velocity` | `HDiv` (vector) | `DIV_filtration_velocity` |

Only **2D and 3D** cell topologies are supported.

### Equation-set parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Model ID` | `""` | Closure model sublist id |
| `Basis Order` | `1` | FE order |
| `Integration Order` | `2` | Quadrature order |
| `Displacement Basis Type` | `HGrad` | Displacement FE space |
| `Pressure Basis Type` | `HGrad` | Pressure FE space |
| `Filtration Basis Type` | `HDiv` | Darcy flux FE space |
| `Inverse Permeability` | `inverse_permeability` | Closure field \(\kappa^{-1}\) |
| `Fluid Source` | `fluid_source` | Closure field \(G_b\) |
| `Body Force` | `body_force` | Closure field \(\mathbf{F}_b\) |
| `Poroelastic Stress Flux` | `poroelastic_stress_flux` | Closure field for solid momentum |
| `Coupled Fluid Traction` | `coupled_fluid_traction` | Interface traction from the fluid side |
| `Structure Density` | `structure_density` | Closure field \(\rho_b\) |
| `Storage Coefficient` | `storage_coefficient` | Closure field \(c_0\) |
| `Biot-Willis Parameter` | `biot_willis` | Closure field \(\alpha\) |
| `Spring Coefficient` | `spring_coefficient` | Closure field \(\gamma\) |

Enable transient terms with `Assembly` → `Build Transient Support: true` in the driver input.

---

## Weak residuals

`EquationSet_Biot::buildAndRegisterEquationSetEvaluators` assembles three residual groups with standard Panzer integrators, then sums sub-operators per DOF.

### Solid momentum (`displacement`)

| Operator | Integrator | Closure / field inputs |
|----------|------------|-------------------------|
| Transient | `Integrator_BasisTimesVector` | `DXDT_displacement` × `structure_density` |
| Pore-pressure coupling | `Integrator_BasisTimesVector` | `GRAD_pore_pressure` × `biot_willis` |
| Poroelastic stress | `Integrator_BasisTimesVector` | `poroelastic_stress_flux` |
| Spring | `Integrator_BasisTimesVector` | `displacement` × `spring_coefficient` |
| Body force | `Integrator_BasisTimesVector` | `body_force` |
| Coupled traction | `Integrator_BasisTimesVector` | `coupled_fluid_traction` |

### Mass conservation (`pore_pressure`)

| Operator | Integrator | Closure / field inputs |
|----------|------------|-------------------------|
| Storage | `Integrator_BasisTimesScalar` | `DXDT_pore_pressure` × `storage_coefficient` |
| Solid dilatation | `Integrator_GradBasisDotVector` | `DXDT_displacement` × `biot_willis` |
| Filtration divergence | `Integrator_BasisTimesScalar` | `DIV_filtration_velocity` |
| Source | `Integrator_BasisTimesScalar` | `fluid_source` |

### Darcy law (`filtration_velocity`)

| Operator | Integrator | Closure / field inputs |
|----------|------------|-------------------------|
| Pressure gradient | `Integrator_DivBasisTimesScalar` | `pore_pressure` |
| Permeability | `Integrator_BasisTimesVector` | `filtration_velocity` × `inverse_permeability` |

---

## Closure model factory

`flujo::ClosureModelFactory` extends `panzer::ClosureModelFactory`. Instantiate it through `ClosureModelFactory_TemplateBuilder` when wiring `Driver` (the driver factory is still a stub in-tree).

### Constant closures

- Scalar: `Value: <double>` → `panzer::Constant` at IP and on each unique basis.
- Vector: `Value X`, `Value Y`, (`Value Z` in 3D) → `panzer::ConstantVector` at IP.

### Typed Biot / Stokes–Biot closures

| `Type` | Evaluator | Purpose |
|--------|-----------|---------|
| `BIOT POROELASTIC STRESS` | `BiotPoroelasticStress` | Fills `poroelastic_stress_flux` from `GRAD_displacement` and `pore_pressure` |
| `STOKES BIOT FLUID TRACTION` | `StokesBiotFluidTraction` | Normal and tangential fluid traction from `pressure` and `GRAD_velocity` |
| `STOKES BIOT NORMAL FILTRATION` | `StokesBiotNormalFiltration` | Scalar jump \((\mathbf{u}\cdot\mathbf{n}) - (\partial_t\eta\cdot\mathbf{n})\) |
| `STOKES BIOT BJS SLIP` | `StokesBiotBeaversJosephSlip` | Tangential slip residual \(\boldsymbol{\tau}_t + \beta \mathbf{u}_t\) |

**`BIOT POROELASTIC STRESS` parameters:** `Lame Mu`, `Lame Lambda`, `Biot-Willis Parameter`, `Displacement Gradient Name` (default `GRAD_displacement`), `Pore Pressure Name` (default `pore_pressure`).

**`STOKES BIOT FLUID TRACTION` parameters:** `Fluid Viscosity`, `Fluid Pressure Name`, `Fluid Velocity Gradient Name`, `Normal Traction Name`, `Tangential Traction Name`, plus interface normal via `Interface Normal` (`X` / `Y` / `Z`) or `Interface Normal 0`, `Interface Normal 1`, … (normalized automatically).

**`STOKES BIOT NORMAL FILTRATION` parameters:** `Fluid Velocity Name`, `Structure Velocity Name` (typically `DXDT_displacement`), and the same interface-normal keys.

**`STOKES BIOT BJS SLIP` parameters:** `Slip Coefficient`, `Fluid Velocity Name`, `Tangential Traction Name`.

The factory also registers non-Biot types (`ARRHENIUS REACTION SOURCE`, `SUPG SCALAR TRANSPORT`) used by other equation sets.

### Poroelastic stress flux note

`BiotPoroelasticStress` evaluates a Saint-Venant–Kirchhoff stress tensor at quadrature points and stores, for each component \(i\), the row sum \(\sum_j \sigma_{ij}\) in `poroelastic_stress_flux`. That vector is paired with `Integrator_BasisTimesVector` as a lumped stress contribution, not a full symmetric-gradient weak form of \(\nabla\cdot\sigma_b\).

---

## Coupling to Navier–Stokes

The **Navier–Stokes** equation set (`Type: Navier-Stokes`) registers `velocity`, `pressure`, and their gradients but does **not** yet assemble fluid residuals. Stokes–Biot coupling is therefore closure-driven: fluid fields are inputs to interface evaluators, and resulting traction or slip fields feed Biot residuals or future fluid BCs.

Typical field wiring on a shared interface:

| Fluid (NS) field | Biot / closure consumer |
|------------------|-------------------------|
| `pressure` | `STOKES BIOT FLUID TRACTION` |
| `GRAD_velocity` | `STOKES BIOT FLUID TRACTION` |
| `velocity` | `STOKES BIOT NORMAL FILTRATION`, `STOKES BIOT BJS SLIP` |
| `DXDT_displacement` | `STOKES BIOT NORMAL FILTRATION` |
| `stokes_biot_normal_traction` (evaluated) | Map into `coupled_fluid_traction` on the solid side |
| `stokes_biot_tangential_traction` | BJS slip closure |

Use the same `Model ID` on both physics blocks when a single closure sublist should build shared interface fields.

### Example closure sublist (YAML)

```yaml
Closure Models:
  biot_stokes_coupling:
    inverse_permeability:
      Value: 1.0e-12
    storage_coefficient:
      Value: 1.0e-4
    biot_willis:
      Value: 1.0
    structure_density:
      Value: 1000.0
    spring_coefficient:
      Value: 0.0
    body_force:
      Value X: 0.0
      Value Y: 0.0
      Value Z: 0.0
    fluid_source:
      Value: 0.0
    poroelastic_stress_flux:
      Type: BIOT POROELASTIC STRESS
      Displacement Gradient Name: GRAD_displacement
      Pore Pressure Name: pore_pressure
      Lame Mu: 1.0e6
      Lame Lambda: 1.0e6
      Biot-Willis Parameter: 1.0
    coupled_fluid_traction:
      Type: STOKES BIOT FLUID TRACTION
      Fluid Pressure Name: pressure
      Fluid Velocity Gradient Name: GRAD_velocity
      Normal Traction Name: stokes_biot_normal_traction
      Tangential Traction Name: stokes_biot_tangential_traction
      Fluid Viscosity: 0.0035
      Interface Normal Z: 1.0
    filtration_jump:
      Type: STOKES BIOT NORMAL FILTRATION
      Fluid Velocity Name: velocity
      Structure Velocity Name: DXDT_displacement
      Interface Normal Z: 1.0
    bjs_slip:
      Type: STOKES BIOT BJS SLIP
      Fluid Velocity Name: velocity
      Tangential Traction Name: stokes_biot_tangential_traction
      Slip Coefficient: 1.0e-4
```

Point both the Biot and Navier–Stokes physics blocks at `Model ID: biot_stokes_coupling` (or split models if you separate material and interface entries).

---

## Related documentation

- [CODE_DOCUMENTATION.md](CODE_DOCUMENTATION.md) — repository-wide architecture, factory dispatch, and implementation status.
- [arXiv:2507.10538](https://arxiv.org/abs/2507.10538) — full Stokes–Biot–plate formulation, splitting scheme, and stability analysis.
