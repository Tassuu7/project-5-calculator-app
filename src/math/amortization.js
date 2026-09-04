/**
 * Financial Amortization & Cash Flow Module
 * Computes monthly loan schedules, principal/interest breakdown, and net present value.
 */

class FinancialAmortization {
  static calculateSchedule(principal, annualRate, years) {
    const monthlyRate = annualRate / 100 / 12;
    const totalMonths = years * 12;
    const monthlyPayment = (principal * monthlyRate * Math.pow(1 + monthlyRate, totalMonths)) /
      (Math.pow(1 + monthlyRate, totalMonths) - 1);

    const schedule = [];
    let balance = principal;

    for (let m = 1; m <= totalMonths; m++) {
      const interest = balance * monthlyRate;
      const principalPaid = monthlyPayment - interest;
      balance = Math.max(0, balance - principalPaid);
      schedule.push({
        month: m,
        payment: monthlyPayment,
        principal: principalPaid,
        interest: interest,
        remainingBalance: balance
      });
    }

    return {
      monthlyPayment,
      totalPayment: monthlyPayment * totalMonths,
      totalInterest: monthlyPayment * totalMonths - principal,
      schedule
    };
  }

  static netPresentValue(rate, cashFlows) {
    return cashFlows.reduce((acc, cf, t) => acc + cf / Math.pow(1 + rate, t), 0);
  }
}

if (typeof module !== 'undefined') module.exports = FinancialAmortization;
