import numpy as np

class Star(object):

  def __init__(self, star_file):
    self._star_file  = star_file
    self._header  = header_reader(star_file)
    self._entries = body_reader(star_file, self._header)

  def get_object(self, objnum, *args):
    self.obj  = objnum
    entry = self._entries[self._entries['NUMBER'] == self.obj]
    if entry.shape[0] != 1:
      print('WARNING!!! More than one entry found for number {}'.\
            format(self.obj))
      print('Please check {}'.format(self._star_file))

    if args:
      return tuple(entry[list(map(str.upper, args))])[0]
    else:
      return tuple(entry)[0]
   
def body_reader(star_file, header):
  body = np.genfromtxt(star_file, names=tuple(header), dtype=None)
  return body

def header_reader(star_file):
  with open(star_file, 'r') as f:
    headernames = [line.strip('#').split()[1] for line in f\
                   if line.startswith('#')]
  header = np.array(headernames)
  return header
