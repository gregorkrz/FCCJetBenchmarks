import argparse

# print the arguments received from command line
import sys
print("Arguments received from command line:", sys.argv)

parser = argparse.ArgumentParser()
parser.add_argument("--test-arg", type=int)

# parse only arguments after the first "--" if present
if "--" in sys.argv:
    argv_after_sep = sys.argv[sys.argv.index("--") + 1:]
    args = parser.parse_args(argv_after_sep)
else:
    args = parser.parse_args()


print("Test arg:", args.test_arg)

def build_graph(df, dataset):
    return [], 0
