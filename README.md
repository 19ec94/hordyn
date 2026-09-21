# Follicular Dynamics
This project studies ovarian follicle development in mammals in response to interactions between various hormones. The complex dynamics between hormones and their effects on follicular development are the central theme of this work.

This repository contains the source code for the simulations presented in the accompanying research paper.

# Bridge between PDE and ODE Models

This repository contains simulations based on PDE models as well as their corresponding ODE models. The application is built using CMake, C++, and Python.

The source code for this project can be found under:

```bash
src/bridge-ode-pde
```

## Requirements

- CMake
- C++
- Python

## How to compile?

From the project root,run:

```bash
mkdir build
cmake -S . -B build
cmake --build build
```
This will configure the project in the `build` directory and then compile the application.

## Run the code

After building the project, run the binary:
```bash
./bin/bridge-pde-ode
```

## Contact

| Name | Email |
| --- | --- |
| Dr. Edilbert Christhuraj | [edilbert.christhuraj@hs-anhalt..de](mailto:edilbert.christhuraj@hs-anhalt.de) |
| Prof. Alexander Lange | [alexander.lange@hs-anhalt.de](mailto:alexander.lange@hs-anhalt.de) |
| Dr. Claudio Iuliano | [claudio.iuliano@hs-anhalt..de](mailto:claudio.iuliano@hs-anhalt.de) |

