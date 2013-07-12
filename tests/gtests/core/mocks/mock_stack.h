class Mock_Stack : public _Stack {
 public:
  MOCK_METHOD0(Initialize,
      void(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
};
