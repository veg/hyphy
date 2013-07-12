class Mock_HBLCommandExtras : public _HBLCommandExtras {
};

class Mock_ElementaryCommand : public _ElementaryCommand {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(toStr,
      BaseRef(void));
};
