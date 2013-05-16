class Mock_AssociativeList : public _AssociativeList {
 public:
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD4(Execute,
      _PMathObj(long, , , _hyExecutionContext));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD0(Compute,
      _PMathObj(void));
  MOCK_METHOD1(Merge,
      void(_PMathObj));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(ObjectClass,
      unsigned long(void));
};
