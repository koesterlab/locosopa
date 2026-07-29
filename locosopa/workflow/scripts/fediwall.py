"""
Snakemake script to download and configure fediwall for the Snakemake community.
"""

import os
import tempfile
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path
from typing import Any, Dict

try:
    import yaml
except ImportError:
    yaml = None


def main():
    # Access Snakemake objects
    name = snakemake.params.name
    hashtags = snakemake.params.hashtags
    accounts = snakemake.params.accounts
    fedi_wall_version = snakemake.params.version
    fedi_wall_sha256 = snakemake.params.expected_sha256
    fediwall_dir = Path(snakemake.output[0])

    # Normalize hashtags: ensure list of lowercase strings without leading '#'
    if not isinstance(hashtags, list):
        raise TypeError("'hashtags' must be a list")
    if not isinstance(accounts, list):
        raise TypeError("'accounts' must be a list")

    # Configuration
    config: Dict[str, Any] = {
        "title": f"Fediwall for {name}",
        "accounts": accounts,
        "hashtags": [str(h).lstrip("#").lower() for h in hashtags],
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
        "showreplies": True,
        "showsensitive": False,
    }

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
            fediwall_dir.mkdir(parents=True, exist_ok=True)

            # Copy fediwall files
            for item in source_dir.iterdir():
                if item.is_dir():
                    shutil.copytree(item, fediwall_dir / item.name, dirs_exist_ok=True)
                else:
                    shutil.copy2(item, fediwall_dir / item.name)

        # Clean up downloaded file
        os.unlink(tmp_file.name)

    # Write merged config file
    with open(fediwall_dir / "wall-config.json", "w") as f:
        json.dump(config, f, indent=2)


if __name__ == "__main__":
    main()
