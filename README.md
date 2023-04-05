# Template Repository

> Folder structure and naming conventions for SERG Group python repositories

### Directory layout:

    .
    ├── calculation  .............................. # Files used for class testing and analysis implementations
    │    ├── 2022-10-04 - Initial Calculation.py .. # Example File
    │    └── ...
    │
    ├── main_code  ................................ # Contain classes and sublcasses 
    │    ├── main_class.py  ....................... # Example File
    │    └── ...
    │
    ├── documentation ............................. # Contain old code and usefull documents  
    │    ├── 2022-10-05 - Test  ................... # Example Folder
    │    └── ...
    │ 
    ├── test  ..................................... # Contain Unittest files
    │    ├── initial_test.py  ..................... # Example File
    │    └── ...
    │
    └── requirements.txt  ......................... # A txt file containing all the libraries needed to run the code

### Naming Conventions:

| Type                             | Naming Convention                                                                                                                                                                              | Examples                                                                 |
|:---------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------|
| **Folder**                       | **lowercase letter** separate words with **underscores**                                                                                                                                       | .../**body_parts**                                                       |
| **File**                         | Same as folder, if contains only one class should have **the same name of the class**                                                                                                          | .../**body_model.py**                                                    |
| **Py Class**                     | Use the python convention: <br/><pre>1. **Start** each word with a **capital letter**.<br/>2. **Do not** separate words with underscores.</pre>                                                | Class **BodyModel**:                                                     |
| **Py Class Method / Properties** | Use the python convention: <br/><pre>1. Use **lower case letter** only.<br/>2. Separate words with **underscores**.<br/>3. For private methods/properties start with a double underscore</pre> | **cylinder_leght**=10<br/>**__cylinder_leght**=10<br/>def **try_this**() |
