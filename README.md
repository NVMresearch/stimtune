# StimTune
An SPM toolbox for optimising the timing of fMRI experiments.

This toolbox helps design fMRI experiments which have the following qualities:

1. A factorial design with any number of factors and levels.
2. Stimuli with the same duration across all trials (which may be an unknown duration).
3. Inter-trial intervals (ITIs) with durations that differ across trials, sampled from a distribution (e.g. uniform, gamma). The parameters of that distribution may be unknown.

The toolbox uses a simple grid search to assess mixtures of stimulus duration and parameters of the ITI sampling distribution. For each mixture of duration and parameters, it generates a series of design matrices with randomised trial order and ITIs. It calculates the statistical efficiency of each design for each contrast, and keeps track of the best designs overall.

See `examples/run_examples.m` to get started.
