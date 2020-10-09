
import os
import zipfile
from github_release import gh_asset_download
import argparse

HERE = os.path.abspath(os.path.dirname(__file__))


def extract(repo="DOE-ICoM/icom-mesh-data", vtag="v0.1.0"):
    """
    EXTRACT: download and unpack data assets for icom-mesh.

    Authors: Darren Engwirda

    """
    cwd_pointer = os.getcwd()

    try:
        print("Extracting to:", os.path.join(HERE, "data"))
        print("")

        directory = os.path.join(HERE, "data")

        os.chdir(directory)

        gh_asset_download(repo, vtag)

        for item in os.listdir(directory):
            if item.endswith(".zip"):
                file_name = os.path.abspath(item)
                zip_func = zipfile.ZipFile(file_name)
                zip_func.extractall(directory)
                zip_func.close()

        print("")
        print("Done downloading + unpacking data assets.")
        print("")

    finally:
        os.chdir(cwd_pointer)

    return


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    REPO_NAME = "DOE-ICoM/icom-mesh-data"

    parser.add_argument(
        "--repo_name", dest="repo_name", type=str,
        required=False,
        help="icom-mesh-data repo. name", default=REPO_NAME)

    parser.add_argument(
        "--repo_vtag", dest="repo_vtag", type=str,
        required=False,
        help="icom-mesh-data release tag", default="v0.1.0")

    args = parser.parse_args()

    extract(repo=args.repo_name, vtag=args.repo_vtag)
