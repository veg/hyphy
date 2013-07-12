class Mock_TreeTopology : public _TreeTopology {
 public:
  MOCK_METHOD1(PreTreeConstructor,
      void(bool));
  MOCK_METHOD2(MainTreeConstructor,
      bool(, ));
  MOCK_METHOD1(PostTreeConstructor,
      void(bool));
  MOCK_METHOD0(FlatRepresentation,
      _PMathObj(void));
  MOCK_METHOD1(toFileStr,
      void(FILE));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(Compare,
      _FString*(_PMathObj));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD6(FinalizeNode,
      bool(node, , , , , ));
  MOCK_METHOD4(Execute,
      _PMathObj(, , , _hyExecutionContext));
  MOCK_METHOD2(EdgeCount,
      void(long &, long));
  MOCK_METHOD0(TipCount,
      _PMathObj(void));
  MOCK_METHOD0(BranchCount,
      _PMathObj(void));
  MOCK_METHOD1(AVLRepresentation,
      _PMathObj(_PMathObj));
  MOCK_METHOD0(ObjectClass,
      unsigned long(void));
  MOCK_METHOD1(FindCOT,
      _AssociativeList*(_PMathObj));
  MOCK_METHOD3(GetNodeName,
      void(node, , ));
  MOCK_METHOD3(GetBranchLength,
      void(node, , ));
  MOCK_METHOD2(GetBranchLength,
      void(node<long> *, _Parameter));
  MOCK_METHOD2(GetBranchValue,
      void(node<long> *, _String));
  MOCK_METHOD3(GetBranchVarValue,
      void(node<long> *, _String &, long));
  MOCK_METHOD4(PasteBranchLength,
      void(node, , , _Parameter));
  MOCK_METHOD1(TipName,
      _PMathObj(_PMathObj));
  MOCK_METHOD3(BranchName,
      _PMathObj(, , ));
  MOCK_METHOD1(BranchLength,
      _PMathObj(_PMathObj));
  MOCK_METHOD1(RerootTree,
      _PMathObj(_PMathObj));
};
