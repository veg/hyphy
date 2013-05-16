class MockFileState : public FileState {
};

class Mock_DSHelper : public _DSHelper {
 public:
};

class Mock_DataSet : public _DataSet {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD3(operator,
      char(unsigned long, unsigned long, unsigned int));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(toFileStr,
      void(FILE *dest));
};
