class Mock_PolynomialData : public _PolynomialData {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseObj*(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
};
