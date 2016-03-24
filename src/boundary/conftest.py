import sys
import pytest

def pytest_addoption(parser):
  parser.addoption("--build_path", action="store", default=None,
                   help='set build directory path'
                  )
  parser.addoption("--N1", action="store", default=64,
                   help='grid zones in X1'
                  )
  parser.addoption("--N2", action="store", default=64,
                   help='grid zones in X2'
                  )
  parser.addoption("--N3", action="store", default=64,
                   help='grid zones in X3'
                  )
  parser.addoption("--dim", action="store", default=64,
                   help='grid dimension'
                  )


def pytest_configure(config):
  buildPath =  config.getvalue('build_path')
  gridPath  = buildPath + '/grid/'
  boundaryPath = buildPath + '/boundary/'
  sys.path.append(gridPath)
  sys.path.append(boundaryPath)
