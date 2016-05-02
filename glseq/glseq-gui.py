import argparse
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs a GlSeq2 bioinformatic pipeline graphic user interface.\n")

    path = __file__.split("/")
    path = path[:len(path)-1]
    string_path = ""
    for p in path:
        string_path += p + "/"

    subprocess.call(["java", "-jar", string_path + "User Interface/GLSeq2_UI.jar"])