# CE1_PSIM
Simulação e análise de circuitos RLC de segunda ordem por meio de modelagem matemática no MATLAB e simulação no PSIM, com comparação entre resultados teóricos e computacionais para diferentes condições de amortecimento.

# RLC Circuit Simulation and Analysis

Modeling, numerical simulation, and comparative analysis of second-order RLC circuits using **MATLAB** and **PSIM**.

## 📖 Overview

This project aims to study the transient response of a second-order RLC circuit through two complementary approaches:

* Mathematical modeling based on the circuit differential equations.
* Numerical simulation using MATLAB.
* Circuit simulation using PSIM.
* Comparison between theoretical and simulated results.

The analysis covers the three damping conditions:

* **Overdamped**
* **Critically Damped**
* **Underdamped**

The project also investigates the influence of the circuit parameters (**R**, **L**, and **C**) and the initial conditions on the system's dynamic behavior.

---

## 🎯 Objectives

* Develop the mathematical model of a second-order RLC circuit.
* Implement the numerical solution in MATLAB.
* Simulate the same circuit in PSIM.
* Compare both approaches.
* Validate the mathematical model through simulation.
* Analyze the transient response for different damping conditions.

---

## 📂 Repository Structure

```text
.
├── MATHLAB/
|   ├── README.md
│   ├── CODE.m
│   ├── VC.png
│   ├── IT.png
│   ├── Criticamente_Amortecido.png
│   ├── Subamortecido.png
│   └── Superamortecido.png
│
├── PSIM/
|   ├── CRITICAMENTE AMORTECIDO
|   |   ├── README.md
|   |   ├── CIRCUITO.png
|   |   ├── IT.png
|   |   ├── VC.png
|   |   ├── VR.png
|   |   ├── VI.png
|   |   ├── MEDIÇÕES.txt
|   |   ├── Criticamente Amortecido.psimsch
│   ├── SUBAMORTECIDO
|   |   ├── README.md
|   |   ├── CIRCUITO.png
|   |   ├── IT.png
|   |   ├── VC.png
|   |   ├── VR.png
|   |   ├── VI.png
|   |   ├── MEDIÇÕES.txt
|   |   ├── Subamortecido.psimsch
│   ├── SUPERAMORTECIDO
|   |   ├── README.md
|   |   ├── CIRCUITO.png
|   |   ├── IT.png
|   |   ├── VC.png
|   |   ├── VR.png
|   |   ├── VI.png
|   |   ├── MEDIÇÕES.txt
|   |   ├── Superamortecido.psimsch
│
├── report/
│   ├── figures/
│   ├── ieee/
│   └── references/
│
├── docs/
│
└── README.md
```

---

## 📊 Simulations

The following cases are analyzed:

* Superdamped (Overdamped)
* Critically Damped
* Underdamped

For each case, the following quantities are evaluated:

* Capacitor voltage
* Inductor current
* Transient response
* Settling time
* Oscillation frequency (when applicable)

---

## 🛠️ Software

* MATLAB
* PSIM
* Overleaf (IEEE Template)

---

## 📄 Report

The accompanying paper follows the IEEE article format and contains the following sections:

* Abstract
* Introduction
* RLC Circuit Analysis
* Results
* Conclusion

All figures and plots are exported in **vector format** for high-quality publication.

---

## 📚 Theory

The project is based on the analysis of second-order linear systems described by:

[
L C \frac{d^2v_c}{dt^2} + R C \frac{dv_c}{dt} + v_c = V_s
]

whose characteristic equation is

[
s^2+\frac{R}{L}s+\frac{1}{LC}=0
]

The damping condition depends on the damping ratio:

* ζ > 1 → Overdamped
* ζ = 1 → Critically Damped
* ζ < 1 → Underdamped

---

## 📈 Expected Results

The numerical solution implemented in MATLAB is expected to closely match the PSIM simulations, validating the mathematical model and illustrating the effect of the circuit parameters on the transient response.

---

## 👨‍💻 Author

**Ricardo Barbosa**

Electrical Engineering Student

---

## 📜 License

This project is intended for academic and educational purposes.
