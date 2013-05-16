class Mock_Site : public _Site {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(Clear,
      void(void));
};
