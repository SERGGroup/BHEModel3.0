from main_code.main_class import Main
import unittest


class MyTestCase(unittest.TestCase):

    def test_something(self):

        text = "trial"
        main = Main(text)
        self.assertEqual(main.text, text)  # add assertion here


if __name__ == '__main__':
    unittest.main()
