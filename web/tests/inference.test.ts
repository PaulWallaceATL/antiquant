import { predictClassical, predictQuantum } from '../lib/inference';

// Mock the spawn function to avoid running actual Python scripts during simple unit tests
// In a real scenario, we might use a dedicated mock library or integration tests.
// For now, we'll just test the TS wrapper logic if possible, or create a simple test file.

describe('Inference API', () => {
    it('should have predictClassical function', () => {
        expect(typeof predictClassical).toBe('function');
    });

    it('should have predictQuantum function', () => {
        expect(typeof predictQuantum).toBe('function');
    });
});
