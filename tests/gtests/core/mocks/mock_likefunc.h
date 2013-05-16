class MockMSTCache : public MSTCache {
};

class Mock_LikelihoodFunction : public _LikelihoodFunction {
 public:
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(Compute,
      _Parameter(void));
  MOCK_METHOD0(Optimize,
      _Matrix*());
  MOCK_METHOD1(CovarianceMatrix,
      _PMathObj());
  MOCK_METHOD0(RescanAllVariables,
      void(void));
  MOCK_METHOD0(ScanAllVariables,
      void(void));
};
