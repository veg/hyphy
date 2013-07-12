class Mock_CalcNode : public _CalcNode {
 public:
  MOCK_METHOD0(Compute,
      _PMathObj(void));
  MOCK_METHOD0(ObjectClass,
      unsigned long(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD1(FreeUpMemory,
      long(long));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD0(HasChanged,
      bool(void));
  MOCK_METHOD1(NeedToExponentiate,
      bool());
  MOCK_METHOD2(SetModel,
      void(long, _AVLListXL));
  MOCK_METHOD1(SetDependance,
      long(long));
  MOCK_METHOD0(RemoveModel,
      void(void));
  MOCK_METHOD2(ReplaceModel,
      void(_String &modelName, _VariableContainer *parentTree));
  MOCK_METHOD0(CheckForReferenceNode,
      long(void));
  MOCK_METHOD0(ClearCategoryMap,
      void(void));
  MOCK_METHOD3(SetupCategoryMap,
      void(_List &, _SimpleList &, _SimpleList));
};

class MocknodeCoord : public nodeCoord {
 public:
};
