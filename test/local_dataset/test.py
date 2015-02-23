class ProgressBar (object):
    def __init__ (self, tot, number_step):
        self.tot = tot
        self.number_step = number_step
        self.numeric_step = int(tot/number_step)
        self.n_step = 0
        

    def __call__ (self, n):
        if n%self.numeric_step == 0:
            if self.n_step == self.number_step:
                print("[{}] 100% DONE".format("X"*self.n_step))
            else:
                print("[{}{}] {}%".format(
                "X"*self.n_step,
                "-"*(self.number_step - self.n_step),
                self.n_step*100/self.number_step))
                
            self.n_step +=1

        
