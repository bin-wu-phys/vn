from multiprocessing import Process, Queue
# import os

class Proc(Process):
    """
    Subclass Process.

    """
    def __init__(self, nrange, rLIST, ob, q):
        Process.__init__(self)
        self.nrange = nrange
        self.rLIST = rLIST
        self.ob = ob
        self.q = q
        # print('Pid of {} in init: {}'.format(self.name, os.getpid()))

    def start(self):
        Process.start(self)
        # print('Pid {} is asigned to {}, which deals with {}.\n'.format(self.pid, self.name, self.nrange))

    def run(self):
        """
        The lesson learn here is that self.res is not really updated.
        """
        for n in self.nrange:
            en, un = self.ob.InterpRestFramen(n, self.rLIST)
            # print('Pid {} has done with {}.\n'.format(self.pid, n))
            self.q.put([n, en, un])
            # print('Pid {} has put the result in the queue.\n'.format(self.pid, n))
        # print('Pid {} has done the job.\n'.format(self.pid))



# class obs:
#     def InterpRestFramen(self, n, rLIST):
#         return n, rLIST*rLIST
#
# if __name__ == '__main__':
#     q = Queue()
#     ps = []
#
#     for n in range(4):
#         p = Proc(n, n, obs(), q)
#         ps.append(p)
#         p.start()
#
#     for p in ps:
#         p.join()
#
#     res = []
#     while not q.empty():
#         ret = q.get()
#         res.append(ret)
#     print(res)
