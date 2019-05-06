class Test:
  def __init__(self,paramIn):
    self.paramIn = paramIn
  def printer(self):
    print self.paramIn
    return self.paramIn[0]
  def other(self):
    print self.paramIn+'haha'
    return self.paramIn+'haha'
