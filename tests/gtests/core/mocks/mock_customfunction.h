class Mock_CustomFunction : public _CustomFunction {
 public:
  MOCK_METHOD0(Compute,
      _Parameter(void));
  MOCK_METHOD0(RescanAllVariables,
      void(void));
};
