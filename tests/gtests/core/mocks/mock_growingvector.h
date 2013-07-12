class Mock_GrowingVector : public _GrowingVector {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(Clear,
      void(void));
  MOCK_METHOD0(GetHDim,
      long(void));
  MOCK_METHOD0(GetVDim,
      long(void));
};
