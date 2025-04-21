# ClosedLoopDBS
 Receive brain recording data and program stimulus in real time.

# Introduction 
 This repository contains programs designed to support closed-loop brain stimulus programmed in real-time using brain recordings using parallel CPU processes. Multiple types of stimulus and recording hardware are supported, including Blackrock and AlphaOmega systems. The repository also contains implementations of two experimental protocols: 
 - Memory: patients with sEEG are tested on their ability to recall pictures, and phase-dependent stimulus is synchronized with the theta brainwave cycle.
 - Parkinson's Disease (PD): patients with cortical (ECoG) and/or deep brain stimulation (DBS) electrodes are directed to plan and execute motor tasks, and phase-dependent stimulus is synchronized with the beta brainwave cycle.

# How to Run 
 ...

# Workflow 
 The following general steps are performed on loop in the background. The first three steps can be executed by iterReadBrain called on loop with appropriate imputs and functions determining the algorithm and protocol. The rest can be executed by ...
 - A batch of recorded raw brain data is obtained from the hardware and organized by channel. Channels can have different sample rates/properties and their data buffers can be different lengths.
 - Incomming data of selected channels is filtered using specified filter coefficients.
 - Data (filtered and/or unfiltered) is modelled and processed; this process may be used to forecast future data, and outputs are calculated that may be used for stimulus (e.g. timing of next pulse).
 - Stimulus is scheduled using any/all of the above data according to a specified controller algorithm.

# Phase-Dependent Stimulation (PDS) 
 The script(s) bg_PhaseDetect can perform PDS as an independent process, accepting parallel pool data queue(s) that report data back to the user and a user data structure that can be the handles or app structure from a MATLAB GUI. The script(s) pollDataQueue_PhaseDetect can be used in an interface to monitor the data queue and organize the data in a plottable format. The script(s) postEval_PDS can report the accuracy of the real-time algorithm offline at a later time. 
 The algorithm monitors a number of channels (e.g. for grid display of power or modulation index) and has one channel designated for recording, which may have an extended buffer. The recording channel will be filtered (e.g. into the beta or theta range), artifact-removed, and a set period of future data will be forecast using a supplied autonomous model (e.g. autoregressive model). At all times, the time until the next phase(s) of interest will be predicted using the instantaneous frequency and phase, and this will be used to schedule stimulus pulses when specified experimental conditions are met. 

# Versions 
 - Version 0 (v0): the first version run successfully in the OR/EMU, built using MATLAB GUIDE; does not use parallel processing, and only supports Blackrock hardware. Due to delays caused by having all processes run on a single loop, timing accuracy and frequency of stimulus actually delivered are limited.
 - Version 1 (v1): updated MATLAB GUIDE-based interfaces that use parallel processing using Blackrock or AlphaOmega hardware.
 - Version 2 (planned) will be rebuilt in the mlapp format.

# Internal/Testing Scripts 
 The following scripts are for internal testing purposes and should not be run generally: 
 - scripts beginning with test_
