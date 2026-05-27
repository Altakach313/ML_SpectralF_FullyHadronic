# Spectral-Function Classifier Model Card

## Included models

This repository includes the Keras models trained with all 21
spectral-function bins:

```text
all_spectral_bins/model_1200.h5
all_spectral_bins/model_1600.h5
all_spectral_bins/model_2000.h5
```

The mass labels denote the gluino signal mass in GeV. Each model is a binary
classifier trained to separate the corresponding gluino-pair signal sample from
the `ttbar` background sample.

## Inputs

The feature order is defined in `feature_schema.json`. The public models use 27
features:

- `SpectralFunc0` through `SpectralFunc20`
- `jetspT0` through `jetspT4`
- `eventHT`

The full intermediate data tables also included additional observables such as
jet masses, missing transverse energy, b-jet multiplicity, `deltaPhi`, `mTbjet`,
and the event-shape `C` parameter, but those are not inputs to the public
all-spectral-function `.h5` models.

## Preprocessing

The training workflow used `StandardScaler`:

1. combine signal and background tables,
2. shuffle with random seed 1,
3. split into train, validation, and test samples,
4. fit the scaler on the training features only,
5. transform validation and test features with that same scaler,
6. train the neural network on the transformed arrays.

The `.h5` files therefore expect standardized inputs with the feature order
given in `feature_schema.json`.

## Architecture

The neural-network architecture used in the training notebooks was:

```text
Dense(64, activation="relu")
Dropout(0.5)
Dense(32, activation="relu")
Dropout(0.5)
Dense(1, activation="sigmoid")
```

The output is a signal-like score between 0 and 1.

## Scope

These models are provided for technical reuse and inspection of the workflow.
They are trained on simulated samples and should be used with the same event
selection, feature definitions, and normalization procedure described in the
README.
