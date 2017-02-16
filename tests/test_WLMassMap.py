"""
Example unit tests for WLMassMap package
"""
import unittest
import desc.wlmassmap

class WLMassMapTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'

    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.wlmassmap.WLMassMap(self.message)
        self.assertEquals(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.wlmassmap.WLMassMap)
        foo = desc.wlmassmap.WLMassMap(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
