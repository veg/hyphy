class Mock_TheTree : public _TheTree {
 public:
  MOCK_METHOD0(HasChanged,
      bool(void));
  MOCK_METHOD0(MarkDone,
      void(void));
  MOCK_METHOD6(FinalizeNode,
      bool(node, , , , , ));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(makeDynamicCopy,
      BaseRef(_String));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD0(ObjectClass,
      unsigned long(void));
  MOCK_METHOD4(Execute,
      _PMathObj(, , , _hyExecutionContext));
  MOCK_METHOD1(TEXTreeString,
      _PMathObj(_PMathObj));
  MOCK_METHOD2(PlainTreeString,
      _PMathObj(_PMathObj, _PMathObj));
  MOCK_METHOD3(GetNodeName,
      void(node, , ));
  MOCK_METHOD3(GetBranchLength,
      void(node, , ));
  MOCK_METHOD2(GetBranchLength,
      void(node<long> *, _Parameter));
  MOCK_METHOD2(GetBranchValue,
      void(node<long> *, _String));
  MOCK_METHOD1(GetBranchSpec,
      _String*(node<long>));
  MOCK_METHOD3(GetBranchVarValue,
      void(node<long> *, _String &, long));
  MOCK_METHOD0(ClearConstraints,
      void(void));
  MOCK_METHOD0(RemoveModel,
      void(void));
  MOCK_METHOD1(PreTreeConstructor,
      void(bool));
  MOCK_METHOD1(PostTreeConstructor,
      void(bool));
};
