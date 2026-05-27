#!/usr/bin/env python3
"""Run inference for a serialized spectral-function classifier.

Input features must already use the same standardization convention as the
training arrays. The public model files do not include a raw-feature scaler.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Any

import pandas as pd


def load_schema(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        schema = json.load(handle)
    return list(schema["default_model_features"])


def load_model(path: Path) -> Any:
    suffix = path.suffix.lower()

    if suffix in {".joblib", ".jl"}:
        import joblib

        return joblib.load(path)

    if suffix in {".pkl", ".pickle"}:
        with path.open("rb") as handle:
            return pickle.load(handle)

    if suffix in {".keras", ".h5"}:
        from tensorflow import keras

        return keras.models.load_model(path)

    if suffix == ".onnx":
        import onnxruntime as ort

        return ort.InferenceSession(str(path))

    raise ValueError(f"Unsupported model format: {path.suffix}")


def predict_scores(model: Any, features: pd.DataFrame) -> Any:
    if model.__class__.__module__.startswith("onnxruntime"):
        input_name = model.get_inputs()[0].name
        outputs = model.run(None, {input_name: features.to_numpy("float32")})
        return outputs[0].reshape(len(features), -1)[:, -1]

    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, -1]

    if hasattr(model, "decision_function"):
        return model.decision_function(features)

    predictions = model.predict(features)
    if getattr(predictions, "ndim", 1) > 1:
        return predictions.reshape(len(features), -1)[:, -1]
    return predictions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, type=Path, help="Serialized model artifact.")
    parser.add_argument("--input", required=True, type=Path, help="Input CSV with model-ready event features.")
    parser.add_argument("--output", required=True, type=Path, help="Output CSV with scores.")
    parser.add_argument(
        "--schema",
        type=Path,
        default=Path(__file__).with_name("feature_schema.json"),
        help="Feature schema JSON.",
    )
    args = parser.parse_args()

    feature_names = load_schema(args.schema)
    data = pd.read_csv(args.input)
    missing = [name for name in feature_names if name not in data.columns]
    if missing:
        raise ValueError(f"Input table is missing required columns: {missing}")

    model = load_model(args.model)
    output = data.copy()
    output["score"] = predict_scores(model, data[feature_names])
    output.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
