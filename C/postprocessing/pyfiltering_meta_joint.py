import numpy as np
import sys

class interval():
    def __init__(self,start,end,pval):
        self.start = start
        self.end = end
        self.pval = pval

    def __lt__(self,temp):
        return self.start < temp.start

    def length(self):
        return self.end - self.start + 1

class cluster():
    def __init__(self,start,end,most_sig_int):
        self.start = start
        self.end = end
        self.most_sig_int = most_sig_int
        self.n_intervals = 1

    def update(self,new_int):
        if new_int.start <= (self.end+1):
            self.end = max(new_int.end,self.end)
            self.n_intervals += 1
            if new_int.pval < self.most_sig_int.pval:
                self.most_sig_int = new_int
            elif new_int.pval == self.most_sig_int.pval and new_int.length() < self.most_sig_int.length():
                self.most_sig_int = new_int
            return 1
        else:
            return 0

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ in "__main__":
    n_intervals = file_len(sys.argv[1]) - 1
    if n_intervals > 0:
        T = np.loadtxt(sys.argv[1],delimiter=',',skiprows=1)
        # Correct shape of T in the edge case of a single row
        if n_intervals==1:
            T = np.reshape(T,(1,T.shape[0]))

        intervals = list()
        for i in xrange(T.shape[0]):
            intervals.append(interval(T[i,1],T[i,1]+T[i,0]-1,T[i,2]))
        intervals.sort()

        clusters = list()
        tmp_cluster = cluster(intervals[0].start,intervals[0].end,intervals[0])
        for i in xrange(1,len(intervals)):
            out = tmp_cluster.update(intervals[i])
            if out==0:
                clusters.append(tmp_cluster)
                tmp_cluster = cluster(intervals[i].start,intervals[i].end,intervals[i])
        clusters.append(tmp_cluster)

        Tc = np.zeros((len(clusters),6))
        for i in xrange(len(clusters)):
            Tc[i,0] = clusters[i].most_sig_int.start
            Tc[i,1] = clusters[i].most_sig_int.end
            Tc[i,2] = clusters[i].most_sig_int.pval
            Tc[i,3] = clusters[i].n_intervals
            Tc[i,4] = clusters[i].start
            Tc[i,5] = clusters[i].end

        np.savetxt(sys.argv[2],Tc,fmt='%d %d %e %d %d %d',
                header='most_sig_int_start_pos,most_sig_int_end_pos,p-val,n_intervals_cluster,cluster_start_pos,cluster_end_pos')
