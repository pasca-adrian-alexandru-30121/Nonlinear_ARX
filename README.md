# 📊 Nonlinear ARX Polynomial Model - Linear Regression Implementation

## 📄 Project Description

This project implements a Nonlinear ARX (AutoRegressive with eXogenous inputs) Polynomial Model using configurable orders na, nb, and polynomial degree m.

The model parameters are estimated using Linear Regression, and the trained model can be used in two ways:

- Prediction Mode: Predicting future outputs given current and past inputs/outputs.
- Simulation Mode: Simulating the system's behavior over time.

## 🔧 Problem Description

The general form of the nonlinear ARX model is:

###  y(k)=f(y(k−1),y(k−2),…,u(k−1),u(k−2),…)+e(k)

Where:

- na = order of the output (number of past outputs used)
- nb = order of the input (number of past inputs used)
- nk = system delay
- m = polynomial degree

## 📐 Structure and Parameters Calculation
The regression matrix (regressors) and parameter vector are built using the Newton Multinomial Expansion, which generates all polynomial combinations of the inputs/outputs up to degree m.

The number of regressor terms is computed based on:
![image](https://github.com/user-attachments/assets/a6a0c102-797b-4b54-a16c-f81c30cfa46e)

## 📊 Model Usage
The trained model can be used in two ways:

- Prediction Mode: One-step ahead prediction using known past data.
- Simulation Mode: Autonomous simulation using only initial conditions and inputs.

## 🚀 Algorithm - Multinomial Term Generation
The multinomial terms are generated using the following recursive process:

1. Initialize COMB_EXPONENTI as an empty list.
2. Set EXP_CURENT vector (size na + nb), with EXP_CURENT[1] = m.
3. Iteratively generate all valid combinations of exponents.
4. Store all combinations in COMB_EXPONENTI.

## ✅ Best Solution
For this specific case study, the best identified configuration was:

- na = 2
- nb = 2
- m = 2
This configuration gave the lowest Mean Squared Error (MSE) on the validation data.
