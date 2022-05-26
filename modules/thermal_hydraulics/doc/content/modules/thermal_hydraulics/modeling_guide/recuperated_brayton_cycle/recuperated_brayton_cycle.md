# Recuperated Brayton Cycle

## Introduction

As shown in the Brayton Cycle modeling guide, and the given input files [closed_brayton_cycle.i](thermal_hydraulics/test/tests/problems/brayton_cycle/closed_brayton_cycle.i) and [open_brayton_cycle.i](thermal_hydraulics/test/tests/problems/brayton_cycle/open_brayton_cycle.i), a Brayton Cycle Power Conversion Unit (PCU) consists of a motor, compressor, turbine, and generator all coupled by a single shaft. Detailed descriptions of the compressor and turbine components used in this example can be found in the previous modeling guide [example](modules/thermal_hydraulics/modeling_guide/brayton_cycle/brayton_cycle.md). In the aforementioned example, a simplified startup transient with a simplified heat source was conducted which demonstrated the Thermal Hydraulics module’s capability to produce torque, power, mass flow rate, and pressure ratios for all PCU components based on the shaft speed.

In this example a different heat source and recuperator were added to the same open Brayton Cycle PCU to demonstrate a more nuanced model along with a more realistic piping structure. Also, a Proportional-Integral-Derivative (PID) controller was added to the motor to conduct a more realistic startup transient of the system.

The main heat source is a cylindrical piece of steel 10 m long with a radius of 15 cm. Internal heat generation is supplied to the heat structure by prescribing a total power of 104 kW. Gas which has been compressed in the compressor, and preheated in the recuperator, passes over the exterior of the heat structure and removes heat via forced convection.

The recuperator is represented by an annular cylindrical heat structure of the same material and width as the main heat source but only 5 m in length and without internal heat generation. Instead, hot exhaust gas leaving the turbine transfers heat to the inside of the recuperator, and cooler gas from the compressor outlet removes heat from the outer surface. This preheats the gas entering the main heat source and improves thermal efficiency of the cycle.  

## Transient Description:

Initially the shaft system is at rest and the working fluid and heat structures are the ambient temperature and pressure. At t = 0 s the motor PID activates and ramps the PCU shaft to 90,000 RPM over the course of 1600 s. When t = 1000 s, heat generation in the main heat structure is activated and linearly increases from 0 kW to a maximum power of 104 kW by t = 8600 s. As the working fluid begins to heat, torque supplied to the shaft from the turbine increases. The motor PID control logic slowly decreases the torque supplied by the motor to 0 N·m once the turbine torque is greater than the torque supplied by the motor. This allows for a graceful transition of PCU control from the motor PID to the working fluid and heat source, and prevents drastically over speeding the system above the rated speed of 96,000 RPM.

## Input File description

The recuperated Brayton Cycle example is executed with the input file [recuperated_brayton_cycle.i](thermal_hydraulics/test/tests/problems/recuperated_brayton_cycle/recuperated_brayton_cycle.i). Shown below
