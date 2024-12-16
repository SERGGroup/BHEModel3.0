import os

CODE_DIR = os.path.dirname(__file__)
BASE_CALCULATION_DIR = os.path.join(os.path.dirname(CODE_DIR), "calculation")
CALCULATION_DIR = os.path.join(BASE_CALCULATION_DIR, "general calculations")
ARTICLE_CALCULATION_DIR = os.path.join(BASE_CALCULATION_DIR, "article calculations")
PROJECT_CALCULATION_DIR = os.path.join(BASE_CALCULATION_DIR, "project calculations")
