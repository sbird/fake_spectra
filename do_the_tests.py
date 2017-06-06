from runtests import Tester
import os.path

tester = Tester(os.path.abspath(__file__), "fake_spectra")

tester.main(sys.argv[1:])
