import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans

import tarfile
import os, sys, argparse

def parse_seq_rmsd_to_dm(st_path_to_file):
  """
    Read a file according to their type

    params:
      st_path_to_file - path to distance file
  """
  ts_lines = []
  if st_path_to_file.endswith(".tgz") or st_path_to_file.endswith(".gz"):
    oj_open_tar = tarfile.open(st_path_to_file, "r:gz")
    for oj_member in oj_open_tar.getmembers():
      oj_open_file = oj_open_tar.extractfile(oj_member)
      if oj_open_file is not None:
        ts_lines    = oj_open_file.read().split("\n")
    oj_open_tar.close()
  else:
    with open(st_path_to_file, "r") as oj_fopened:
      if oj_fopened is not None:
        ts_lines    = oj_fopened.read().split("\n")
    oj_fopened.close()
  return read_content(ts_lines)

def read_content(ts_lines = []):
  """
    make distance matrix from precomputed distances

    params:
      ts_lines: list of line
  """
  in_row_count  = 0
  ts_idxes      = []
  in_rows     = int(np.sqrt(len(ts_lines)-1))
  ar_dm = np.zeros([in_rows, in_rows])
  print("INFO: constructed a matrix of dimension {}x{} from {} rows.".format(in_rows, in_rows, len(ts_lines)-1))

  st_exp_no   = "" 
  idx, jdx = 0, 0
  for st_line in ts_lines:
    if not st_line: continue
    ts_columns    = st_line.split(" ")
    in_num_column = len(ts_columns)
    st_pair1      = ts_columns[0] + "_" + ts_columns[1]
    if in_row_count == 3*in_rows: in_row_count = 0 

    st_pair2 = ""
    if in_num_column == 6:
      if in_row_count < 2*in_rows and in_row_count >= in_rows:  
        st_pair2 = "E02_" + ts_columns[2]
      elif in_row_count < 3 * in_rows and in_row_count >= 2*in_rows:  
        st_pair2 = "E03_" + ts_columns[2]
      else:
        st_pair2 = "E01_" + ts_columns[2]
    elif in_num_column == 7:
      st_pair2 = ts_columns[2] + "_" + ts_columns[3]
    else:
      print("Invalid input file...")
      sys.exit(0)

    if not st_pair1 in ts_idxes: ts_idxes.append(st_pair1)
    idx = ts_idxes.index(st_pair1) 

    if not st_pair2 in ts_idxes: ts_idxes.append(st_pair2)
    jdx = ts_idxes.index(st_pair2) 

    # fill distance matrix
    if in_num_column == 6:
      ar_dm[idx, jdx] = float(ts_columns[5])
    elif in_num_column == 7:
      ar_dm[idx, jdx] = float(ts_columns[6])
    in_row_count += 1

  print("Distance matrix size: {}".format(ar_dm.shape))
  return ar_dm
  
def cluster_dm(ar_dm, fo_rmsd, oj_opened):
  """
    Cluster given distance matrix using predefined threshold

    params:
      ar_dm:    distance matrix (nxn) 
      fo_rmsd:  threshold for cluster
  """
  oj_cluster = DBSCAN(eps=fo_rmsd, metric="precomputed").fit(ar_dm)
  ts_unique, ts_count = np.unique(np.array(oj_cluster.labels_), return_counts=True)

  # Write cluster members
  for cno, ccount in zip(ts_unique, ts_count):
    print("# of cluster {}: {} (radius {})".format(cno, ccount,fo_rmsd))
    oj_opened.write(("Clustering radius:{}\n").format(str(fo_rmsd)))
    for in_clust_num in ts_unique:
      ts_members  = [i for i, e in enumerate(oj_cluster.labels_) if e == in_clust_num]
      st_member   = ",".join(map(str, ts_members))
      oj_opened.write(("{}:{}\n").format(str(in_clust_num), st_member))
    oj_opened.write("---\n")

if __name__ == "__main__":
  parser = argparse.ArgumentParser(prog='SegmentCluster.py', description='Segment clustering')
  parser.add_argument('-i','--input_file', required=True, help='Directory path to file')
  oj_args = parser.parse_args(sys.argv[1:])

  ar_dm     = parse_seq_rmsd_to_dm(oj_args.input_file)
  ts_rmsds  = np.arange(1.0, 2.25, 0.25)

  if os.path.exists("ClusterMembers.txt"): os.remove("ClusterMembers.txt")

  with open("ClusterMembers.txt", "w") as oj_opened:
    for fo_rmsd in ts_rmsds: cluster_dm(ar_dm, fo_rmsd, oj_opened)
