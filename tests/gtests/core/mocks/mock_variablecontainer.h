class Mock_VariableContainer : public _VariableContainer {
 public:
  MOCK_METHOD2(SetModel,
      void(long, _AVLListXL));
  MOCK_METHOD0(MarkDone,
      void(void));
  MOCK_METHOD0(IsContainer,
      bool(void));
  MOCK_METHOD0(HasChanged,
      bool(void));
  MOCK_METHOD1(NeedToExponentiate,
      bool());
  MOCK_METHOD4(ScanForVariables,
      void(, , _AVLListX, long));
  MOCK_METHOD2(ScanForDVariables,
      void(_AVLList &, _AVLList));
  MOCK_METHOD4(ScanForGVariables,
      void(, , _AVLListX, long));
  MOCK_METHOD1(IsModelVar,
      bool(long));
  MOCK_METHOD0(IsConstant,
      bool(void));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(RemoveDependance,
      bool(long));
  MOCK_METHOD1(SetDependance,
      long(long));
  MOCK_METHOD0(ClearConstraints,
      void(void));
  MOCK_METHOD1(GetIthIndependent,
      _Variable*(long));
  MOCK_METHOD1(GetIthDependent,
      _Variable*(long));
  MOCK_METHOD1(GetIthParameter,
      _Variable*(long));
  MOCK_METHOD1(CompileListOfDependents,
      void(_SimpleList));
};
