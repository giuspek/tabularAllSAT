(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A3))) (let ((.def_1 (and A3 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A11))) (let ((.def_4 (or A22 .def_3))) (let ((.def_5 (and .def_4 .def_2))) (let ((.def_6 (and A15 A0))) (let ((.def_7 (not A26))) (let ((.def_8 (or A1 .def_7))) (let ((.def_9 (and .def_8 .def_6))) (let ((.def_10 (or .def_9 .def_5))) (let ((.def_11 (not A9))) (let ((.def_12 (not A23))) (let ((.def_13 (or .def_12 .def_11))) (let ((.def_14 (not .def_13))) (let ((.def_15 (not A8))) (let ((.def_16 (or .def_15 .def_12))) (let ((.def_17 (not .def_16))) (let ((.def_18 (= .def_17 .def_14))) (let ((.def_19 (not .def_18))) (let ((.def_20 (or A24 A15))) (let ((.def_21 (= A14 A6))) (let ((.def_22 (not .def_21))) (let ((.def_23 (and .def_22 .def_20))) (let ((.def_24 (not .def_23))) (let ((.def_25 (and .def_24 .def_19))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and .def_26 .def_10))) (let ((.def_28 (not .def_27))) (let ((.def_29 (not A20))) (let ((.def_30 (not A7))) (let ((.def_31 (and .def_30 .def_29))) (let ((.def_32 (not .def_31))) (let ((.def_33 (not A29))) (let ((.def_34 (not A5))) (let ((.def_35 (= .def_34 .def_33))) (let ((.def_36 (not .def_35))) (let ((.def_37 (or .def_36 .def_32))) (let ((.def_38 (not .def_37))) (let ((.def_39 (or A21 A12))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A18))) (let ((.def_42 (or .def_41 A12))) (let ((.def_43 (not .def_42))) (let ((.def_44 (or .def_43 .def_40))) (let ((.def_45 (or .def_44 .def_38))) (let ((.def_46 (or .def_7 A13))) (let ((.def_47 (not .def_46))) (let ((.def_48 (not A25))) (let ((.def_49 (and .def_48 A17))) (let ((.def_50 (and .def_49 .def_47))) (let ((.def_51 (and A2 A6))) (let ((.def_52 (not .def_51))) (let ((.def_53 (or A5 A0))) (let ((.def_54 (not .def_53))) (let ((.def_55 (and .def_54 .def_52))) (let ((.def_56 (and .def_55 .def_50))) (let ((.def_57 (and .def_56 .def_45))) (let ((.def_58 (and .def_57 .def_28))) (let ((.def_59 (not .def_58))) (let ((.def_60 (not A15))) (let ((.def_61 (or .def_60 A1))) (let ((.def_62 (and A9 .def_34))) (let ((.def_63 (= .def_62 .def_61))) (let ((.def_64 (not .def_63))) (let ((.def_65 (and .def_48 A3))) (let ((.def_66 (not A10))) (let ((.def_67 (and .def_66 A12))) (let ((.def_68 (not .def_67))) (let ((.def_69 (or .def_68 .def_65))) (let ((.def_70 (or .def_69 .def_64))) (let ((.def_71 (not .def_70))) (let ((.def_72 (or A3 .def_0))) (let ((.def_73 (and A23 A24))) (let ((.def_74 (or .def_73 .def_72))) (let ((.def_75 (not A4))) (let ((.def_76 (or .def_75 A1))) (let ((.def_77 (and A21 .def_60))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or .def_78 .def_76))) (let ((.def_80 (= .def_79 .def_74))) (let ((.def_81 (and .def_80 .def_71))) (let ((.def_82 (not .def_81))) (let ((.def_83 (not A0))) (let ((.def_84 (and .def_83 A2))) (let ((.def_85 (or .def_33 .def_7))) (let ((.def_86 (or .def_85 .def_84))) (let ((.def_87 (not A6))) (let ((.def_88 (or .def_0 .def_87))) (let ((.def_89 (not A14))) (let ((.def_90 (not A2))) (let ((.def_91 (and .def_90 .def_89))) (let ((.def_92 (or .def_91 .def_88))) (let ((.def_93 (and .def_92 .def_86))) (let ((.def_94 (or A21 .def_11))) (let ((.def_95 (not A13))) (let ((.def_96 (or .def_95 A9))) (let ((.def_97 (and .def_96 .def_94))) (let ((.def_98 (not .def_97))) (let ((.def_99 (and A17 A24))) (let ((.def_100 (not .def_99))) (let ((.def_101 (and A20 .def_0))) (let ((.def_102 (not .def_101))) (let ((.def_103 (and .def_102 .def_100))) (let ((.def_104 (not .def_103))) (let ((.def_105 (or .def_104 .def_98))) (let ((.def_106 (not .def_105))) (let ((.def_107 (or .def_106 .def_93))) (let ((.def_108 (or .def_107 .def_82))) (let ((.def_109 (not .def_108))) (let ((.def_110 (or .def_109 .def_59))) .def_110))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
