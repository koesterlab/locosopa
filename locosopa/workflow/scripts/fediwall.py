"""
Snakemake script to download and configure fediwall for the Snakemake community.
"""

import tempfile
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path
from typing import Any, Dict

try:
    import yaml  # Optional; only needed if user supplies YAML config
except ImportError:  # pragma: no cover - safe fallback
    yaml = None


def main():
    # Access Snakemake objects
    output_config = snakemake.output.config
    fedi_wall_version = snakemake.params.version
    fedi_wall_sha256 = snakemake.params.expected_sha256
    config_source_path = Path(snakemake.params.config_source)

    if not config_source_path.exists():
        raise FileNotFoundError(f"Fediwall config file not found: {config_source_path}")

    # Load user-provided wall config (JSON or YAML). Must contain title, subtitle, hashtags.
    wall_config_user: Dict[str, Any]
    if config_source_path.suffix.lower() in {".yaml", ".yml"}:
        if yaml is None:
            raise RuntimeError(
                "PyYAML not installed but a YAML config was provided. Install pyyaml or use JSON."
            )
        wall_config_user = yaml.safe_load(config_source_path.read_text()) or {}
    else:
        wall_config_user = json.loads(config_source_path.read_text())

    # Basic validation / defaults
    required_top_level = ["title", "subtitle", "hashtags"]
    missing = [k for k in required_top_level if k not in wall_config_user]
    if missing:
        raise ValueError(
            f"Missing required keys in fediwall config: {missing}. Required: {required_top_level}"
        )

    # Provide defaults for optional keys if absent
    defaults: Dict[str, Any] = {
        "accounts": [],
        "servers": [
            "mastodon.social",
            "fosstodon.org",
            "fediscience.org",
            "scholar.social",
            "genomic.social",
        ],
        "theme": "auto",
        "autorefresh": 30,
        "showboosts": True,
        "showreplies": False,
        "showsensitive": False,
    }
    # Merge defaults without overwriting user-provided values
    for k, v in defaults.items():
        wall_config_user.setdefault(k, v)

    # Normalize hashtags: ensure list of lowercase strings without leading '#'
    hashtags = wall_config_user["hashtags"]
    if not isinstance(hashtags, list):
        raise TypeError("'hashtags' must be a list")
    wall_config_user["hashtags"] = [str(h).lstrip("#").lower() for h in hashtags]

    # Download fediwall with checksum verification
    download_url = f"https://github.com/defnull/fediwall/releases/download/{fedi_wall_version}/fediwall-1.4.0.tgz"

    with tempfile.NamedTemporaryFile(suffix=".tgz", delete=False) as tmp_file:
        urllib.request.urlretrieve(download_url, tmp_file.name)

        # Verify SHA256 checksum
        with open(tmp_file.name, "rb") as f:
            actual_sha256 = hashlib.sha256(f.read()).hexdigest()

        if actual_sha256 != fedi_wall_sha256:
            raise ValueError(
                f"SHA256 checksum verification failed! Expected: {fedi_wall_sha256}, Actual: {actual_sha256}"
            )

        print("SHA256 checksum verification passed")

        # Extract to temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            with tarfile.open(tmp_file.name, "r:gz") as tar:
                tar.extractall(temp_dir)

            # Find the extracted directory (strip-components equivalent)
            extracted_dirs = [d for d in Path(temp_dir).iterdir() if d.is_dir()]
            if len(extracted_dirs) == 1:
                source_dir = extracted_dirs[0]
            else:
                # If multiple dirs or files, assume the temp_dir is the source
                source_dir = Path(temp_dir)

            # Create fediwall directory in build
            fediwall_dir = Path("build/fediwall")
            fediwall_dir.mkdir(parents=True, exist_ok=True)

            # Copy fediwall files
            for item in source_dir.iterdir():
                if item.is_dir():
                    shutil.copytree(item, fediwall_dir / item.name, dirs_exist_ok=True)
                else:
                    shutil.copy2(item, fediwall_dir / item.name)

        # Clean up downloaded file
        import os

        os.unlink(tmp_file.name)

    # Write merged config file
    with open(output_config, "w") as f:
        json.dump(wall_config_user, f, indent=2)


if __name__ == "__main__":
    main()
