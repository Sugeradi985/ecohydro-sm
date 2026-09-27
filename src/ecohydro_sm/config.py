"""配置加载：YAML -> 层级字典。

配置里**不放任何密钥**。ASF / Earthdata / CDS 的凭据走环境变量或
~/.netrc，见 docs/workflow.md。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

__all__ = ["load_config", "save_config", "default_config_path"]


def default_config_path() -> Path:
    return Path(__file__).resolve().parents[2] / "configs" / "default.yaml"


def load_config(path: str | Path) -> dict[str, Any]:
    """读取 YAML 配置。"""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"配置文件不存在：{path}")
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    return dict(cfg)


def save_config(cfg: dict, path: str | Path) -> None:
    """写入 YAML 配置。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
