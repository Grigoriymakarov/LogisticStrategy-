# Chemical Transport Network Optimization

## Overview

This project implements a Mixed Integer Linear Programming (MILP) model for the optimization of a chemical transportation network over a five-year planning horizon.

The model determines the optimal fleet composition, vehicle allocation, route utilization, and delivery strategy while minimizing total operational costs and satisfying customer demand.

The optimization is formulated in Python using the PuLP library and solved using the CBC solver.

---

## Problem Description

A logistics company distributes two chemical products:

* **Acid**
* **Base**

Deliveries are performed from a depot located in Liège to multiple customer destinations.

The model considers:

* Customer demand for acid and base products
* Vehicle acquisition and disposal decisions
* Fleet allocation between transport configurations
* Route selection and scheduling
* Vehicle capacity limitations
* Driver working-time constraints
* Vehicle reconfiguration operations
* Fuel consumption costs

The objective is to identify the least-cost distribution strategy over a 5-year planning horizon.

---

## Network Structure

### Depot

* Liège (Belgium)

### Customer Destinations

| Destination | Description |
| ----------- | ----------- |
| A           | Customer A  |
| B           | Customer B  |
| C           | Customer C  |
| G           | Customer G  |
| H           | Customer H  |

The model uses predefined distances:

* Depot → Customer
* Customer → Customer

to generate feasible delivery routes.

---

## Model Features

### Route Enumeration

All feasible routes containing up to three customer stops are generated automatically.

For each route, the model computes:

* Total travel distance
* Travel time
* Optimal customer visiting order

### Fleet Management

Two truck categories are considered:

* Type 1 trucks
* Type 2 trucks

The optimization decides:

* Number of vehicles purchased
* Number of vehicles sold
* Vehicle allocation to different transport configurations

### Product Configurations

Vehicles can operate in multiple configurations:

* Acid transport
* Base transport
* Mixed transport

Vehicle reconfiguration requires downtime and is explicitly modeled.

### Demand Satisfaction

The model guarantees that yearly demand is satisfied for:

* Acid deliveries to each customer
* Base deliveries

---

## Objective Function

The objective minimizes:

* Vehicle acquisition costs
* Fleet operating costs
* Fuel costs
* Reconfiguration penalties

while accounting for vehicle resale value.

---

## Mathematical Formulation

The optimization includes constraints for:

1. Fleet balance
2. Fleet partitioning
3. Vehicle reconfiguration
4. Acid demand satisfaction
5. Base demand satisfaction
6. Vehicle capacity limits
7. Minimum delivery quantities
8. Driver working-time limitations

The resulting model is a Mixed Integer Linear Program (MILP).

---

## Technologies Used

* Python
* PuLP
* CBC Solver

---

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/chemical-transport-optimization.git
cd chemical-transport-optimization
```

Install dependencies:

```bash
pip install pulp
```

---

## Running the Model

Execute:

```bash
python main.py
```

The program will:

1. Generate all feasible routes
2. Build the optimization model
3. Solve the MILP problem
4. Display:

   * Fleet composition by year
   * Active routes
   * Delivered quantities
   * Total cost
   * Optimization status
   * Demand coverage
   * Vehicle working-time utilization

---

## Note

Developed as part of an operations research and logistics optimization project using Mixed Integer Linear Programming techniques.
