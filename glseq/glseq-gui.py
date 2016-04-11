import argparse
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs a GlSeq2 bioinformatic pipeline graphic user interface.\n")

    subprocess.call(["java", "-jar","User Interface/GLSeq2_UI.jar"])